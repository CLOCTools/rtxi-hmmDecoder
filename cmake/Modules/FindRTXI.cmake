# - Find RTXI
# RealTime eXperimental Interface

# copied from /usr/local/share/rtxi/Makefile.plugin_compile:
set(RTXI_CXX_FLAGS "-pipe -DQT_XML_LIB -DQT_GUI_LIB -DQT_CORE_LIB -DQT_WIDGETS_LIB -DQT_SVG_LIB -DQT_SHARED -Wall -Wextra -g -std=c++11 -fPIC")
set(RTXI_INCLUDE_DIRS "/usr/local/include/rtxi /usr/local/include/rtxi/plugins /usr/local/include/rtxi/libs /usr/local/include /usr/include/x86_64-linux-gnu/qt5/QtOpenGL /usr/include/x86_64-linux-gnu/qt5 /usr/include/x86_64-linux-gnu/qt5/QtPrintSupport /usr/include/x86_64-linux-gnu/qt5 /usr/include/x86_64-linux-gnu/qt5/QtXml /usr/include/x86_64-linux-gnu/qt5 /usr/include/x86_64-linux-gnu/qt5/QtSvg /usr/include/x86_64-linux-gnu/qt5 /usr/include/x86_64-linux-gnu/qt5/QtWidgets /usr/include/x86_64-linux-gnu/qt5 /usr/include/x86_64-linux-gnu/qt5/QtGui /usr/include/x86_64-linux-gnu/qt5 /usr/include/x86_64-linux-gnu/qt5/QtNetwork /usr/include/x86_64-linux-gnu/qt5 /usr/include/x86_64-linux-gnu/qt5/QtCore /usr/include/x86_64-linux-gnu/qt5 /usr/X11R6/include /usr/include/hdf5/serial /usr/include/qwt /usr/local/include/qwt")
set(RTXI_LIBRARIES "Qt5OpenGL Qt5PrintSupport Qt5Xml Qt5Svg Qt5Widgets Qt5Gui Qt5Network Qt5Core git2 qwt-qt5 hdf5 hdf5_hl GL pthread gsl gslcblas m dl alchemy copperplate /usr/xenomai/lib/xenomai/bootstrap-pic.o cobalt modechk pthread rt")
set(RTXI_LIBRARY_DIRS "/usr/local/lib/rtxi/libs /usr/local/lib /usr/lib/x86_64-linux-gnu/hdf5/serial /usr/xenomai/lib")
set(RTXI_SHARED_LINKER_FLAGS "-Wl,--no-as-needed -Wl,@/usr/xenomai/lib/modechk.wrappers")

# cmake needs semi-colon separated lists:
# string(REPLACE " " ";" RTXI_CXX_FLAGS "${RTXI_CXX_FLAGS}")
string(REPLACE " " ";" RTXI_INCLUDE_DIRS "${RTXI_INCLUDE_DIRS}")
string(REPLACE " " ";" RTXI_LIBRARIES "${RTXI_LIBRARIES}")
string(REPLACE " " ";" RTXI_LIBRARY_DIRS "${RTXI_LIBRARY_DIRS}")
# string(REPLACE " " ";" RTXI_SHARED_LINKER_FLAGS "${RTXI_SHARED_LINKER_FLAGS}")

# requires Qt 
find_package(Qt5 REQUIRED COMPONENTS Core)
