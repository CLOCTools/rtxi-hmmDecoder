/*
 * Copyright (C) 2011 Georgia Institute of Technology, University of Utah,
 * Weill Cornell Medical College
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * This is a template implementation file for a user module derived from
 * DefaultGUIModel with a custom GUI.
 */

//printing vector options https://stackoverflow.com/questions/10750057/how-to-print-out-the-contents-of-a-vector

//vector printing code
/*
  printf("\nspikes:\n");
  for (auto i: spike_buff) { printf("%d`",i); }
  printf("\n---\n");

  for (int i=0; i<bufflen; i++)
  {
       printf("%d,",guessed[i]);
  }
  printf("\ndecode done\n");
*/

#include "hmm_decoder.h"
#include "moc_hmm_decoder.cpp"

#include <iostream>
#include <main_window.h>

extern "C" Plugin::Object *
createRTXIPlugin(void)
{
  return new HmmDecoder();
}

static DefaultGUIModel::variable_t vars[] = {
    {
        "spike in",
        "?",
        DefaultGUIModel::INPUT,
    },
    {
        "state out",
        "?",
        DefaultGUIModel::OUTPUT,
    },
    {
        "spikes in>out",
        "?",
        DefaultGUIModel::OUTPUT,
    },

    {
        "FR 1",
        "Firing rate",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "FR 2",
        "Firing rate",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "TR 1",
        "Transition rate",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "TR 2",
        "Transition rate",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "A State",
        "delete me",
        DefaultGUIModel::STATE,
    },
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

HmmDecoder::HmmDecoder(void)
    : DefaultGUIModel("HmmDecoder with Custom GUI", ::vars, ::num_vars), spike_current(0.0), doSample(false)
{
  setWhatsThis("<p><b>HmmDecoder:</b><br>QWhatsThis description.</p>");
  DefaultGUIModel::createGUI(vars,
                             num_vars); // this is required to create the GUI
  customizeGUI();

  initParameters();
  update(INIT); // this is optional, you may place initialization code directly
                // into the constructor
  refresh();    // this is required to update the GUI with parameter and state
                // values
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

HmmDecoder::~HmmDecoder(void)
{
}

void HmmDecoder::execute(void)
{
  // TODO: temporary fudge.
  if (doSample)
  {
    spike_current = ceil(2 * input(0));
  }
  else
  {
    spike_current = 0.0;
  }
  doSample = !doSample;

  //pull from input(0) into buffer
  //decode HMM state in existing buffer
  advanceSpkBuffer(spike_current);
  decodeSpkBuffer();

  //very convoluted. must fix!
  output(0) = state_guess_buff.back();
  output(1) = spike_current;
  //output(0) = ((state_guess_buff.back()<1) ? 0 : 1); //candidate for decoder lag issue
  //-1 is for 0 indexing convention

  return;
}

void HmmDecoder::buildBigHMM()
{
  vFr = {pfr1, pfr2};
  vTr = {ptr1, ptr2};
  trs = vTr;
  frs = vFr;

  /*
    double ptr1_ = (1.0-(ptr1*(nstates-1)));
    double ptr2_ = (1.0-(ptr2*(nstates-1)));

    double pfr1_ = (1.0-pfr1)/(nevents-1);
    double pfr2_ = (1.0-pfr2)/(nevents-1);

    trs = {{ptr1_, ptr1,ptr1}, {ptr1,ptr1_,ptr1}, {ptr1,ptr1,ptr1_}};
    frs = {{pfr1,pfr1_,pfr1_}, {pfr2_,pfr2,pfr2_}, {pfr1_,pfr1_,pfr1}};
*/
  //   =  {             .9  }
  // trs = {{ptr1_, ptr1},{ptr1,ptr1_}};
  //frs = {{20,1,1}, {1,1,20}};
}

void HmmDecoder::initParameters(void)
{

  some_parameter = 0;
  some_state = 0;

  /*
  pfr1=10e-3;
  pfr2=30e-3;
  ptr1=2e-4;
  ptr2=2e-4;
  */
  nstates = 2;
  nevents = 2;
  pfr1 = 1e-3;  //1-1e-2;//
  pfr2 = 20e-3; //.7;//

  ptr1 = 4e-4;
  ptr2 = 4e-4;

  buffi = 0;
  bufflen = 1e3; //300;#default///  3000 is dangerous is tdt room...//

  // [BugFixed] I was tempted to use vector initialization code here, but it was overriding the scope of the vector!
  //vFr.resize(2,0);
  //vTr.resize(2,0);

  spike_buff.resize(bufflen, 1);
  state_guess_buff.resize(bufflen, 0); //unnecessary?

  vFr = {pfr1, pfr2};
  vTr = {ptr1, ptr2};
  //printf("\ngogogogo\n");
  buildBigHMM();
  restartHMM();
  decodeSpkBuffer();

  std::cout << "Decoder!:";

  guess_hmm.printMyParams();
  std::cout << "bufflen{" << bufflen << "}";
}

void HmmDecoder::advanceSpkBuffer(int newSpk)
{

  if (newSpk == 99)
  {
    spike_buff[bufflen - 1] = 99;
    printStuff();
  }
  else
  {
    //cycle buffer left: http://en.cppreference.com/w/cpp/algorithm/rotate
    std::rotate(spike_buff.begin(), spike_buff.begin() + 1, spike_buff.end());
    spike_buff[bufflen - 1] = newSpk;
    //spike_buff[bufflen]=newSpk; // experimental, feel free to delete
    //spike_buff.push(newSpk); //adds to the end
    //spike_buff.pop();
  }
}

int *HmmDecoder::decodeHMM(HMMv guess_hmm_)
{
  int *guessed = viterbi(guess_hmm_, spike_buff, bufflen);
  return guessed;
}

void HmmDecoder::decodeSpkBuffer()
{
  int *guessed = decodeHMM(guess_hmm); //sufficient to cause freeze
  //NB: no idea why this temporary vector is necessary. should be able to replace this with one line...
  std::vector<int> temp_vec(guessed, guessed + bufflen);
  state_guess_buff = temp_vec;

  delete[] guessed; //this closes seq which is dynamically allocated inside viterbi
}

void HmmDecoder::restartHMM()
{
  //really,and internalize parameter modifications from GUI
  //do I actually want to reset the spike buffer? probably not?
  std::vector<double> PI(nstates, .5);

  guess_hmm = HMMv(nstates, nevents, trs, frs, PI);
  //decodeSpkBuffer();//?
}

void HmmDecoder::update(DefaultGUIModel::update_flags_t flag)
{
  int nnn;
  switch (flag)
  {
  case INIT:
    period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
    period_ms = period * 1e-3;
    setState("A State", some_state);

    setParameter("FR 1", pfr1 / period_ms);
    setParameter("FR 2", pfr2 / period_ms);
    setParameter("TR 1", ptr1 / period_ms);
    setParameter("TR 2", ptr2 / period_ms);

    break;

  case MODIFY:
    some_parameter = getParameter("GUI label").toDouble();

    //Need to add the *period*1e3 in here;
    pfr1 = getParameter("FR 1").toDouble() * period_ms;
    pfr2 = getParameter("FR 2").toDouble() * period_ms;
    ptr1 = getParameter("TR 1").toDouble() * period_ms;
    ptr2 = getParameter("TR 2").toDouble() * period_ms;

    vFr = {pfr1, pfr2};
    vTr = {ptr1, ptr2};

    restartHMM();

    //decodeSpkBuffer();
    break;

  case UNPAUSE:
    break;

  case PAUSE:
    break;

  case PERIOD:
    period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
    break;

  default:
    break;
  }
}

void HmmDecoder::customizeGUI(void)
{
  QGridLayout *customlayout = DefaultGUIModel::getLayout();

  QGroupBox *button_group = new QGroupBox;

  QPushButton *abutton = new QPushButton("Button A"); //todo deleteme
  QPushButton *bbutton = new QPushButton("Button B"); //todo deleteme
  QHBoxLayout *button_layout = new QHBoxLayout;
  button_group->setLayout(button_layout);
  button_layout->addWidget(abutton);                                       //DEL
  button_layout->addWidget(bbutton);                                       //DEL
  QObject::connect(abutton, SIGNAL(clicked()), this, SLOT(aBttn_event())); //DEL
  QObject::connect(bbutton, SIGNAL(clicked()), this, SLOT(bBttn_event())); //DEL

  customlayout->addWidget(button_group, 0, 0);
  setLayout(customlayout);
}

void HmmDecoder::printStuff(void)
{
  printf("\n\n decoder;spike_buff=[");
  for (int i = 0; i < bufflen; i++)
  {
    printf("%d,", spike_buff[i]);
  }

  printf("];\nstate_guess=[");
  for (int i = 0; i < bufflen; i++)
  {
    printf("%d,", state_guess_buff[i]);
  }
  printf("];enddec;");
}

// functions designated as Qt slots are implemented as regular C++ functions
void HmmDecoder::aBttn_event(void)
{
  printStuff();
}

void HmmDecoder::bBttn_event(void)
{
  printf(",,");
}
