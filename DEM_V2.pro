TEMPLATE = app
CONFIG += console

CONFIG += opengl
LIBS+=-framework OpenGL -framework GLUT #-ltiff -I/usr/local/include/tiff

# Don't use the following as suggested by Qt when making a new project
#QT -= gui
#TEMPLATE = app
#CONFIG += c++11 console
#CONFIG -= app_bundle
# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
#DEFINES += QT_DEPRECATED_WARNINGS
# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    src/box.cpp \
    src/cell.cpp \
    src/config.cpp \
    src/contact.cpp \
    src/event.cpp \
    src/grain.cpp \
    src/grid.cpp \
    src/in_out.cpp \
    src/main.cpp \
    src/mesh.cpp \
    src/node.cpp \
    src/parameter.cpp \
    src/profile.cpp \
    src/run.cpp \
    src/run_post_process.cpp \
    src/run_visu.cpp \
    src/slice.cpp \
    src/vector.cpp \
    voro/c_loops.cc \
    voro/cell.cc \
    voro/common.cc \
    voro/container.cc \
    voro/container_prd.cc \
    voro/pre_container.cc \
    voro/unitcell.cc \
    voro/v_base.cc \
    voro/v_base_wl.cc \
    voro/v_compute.cc \
    voro/voro++.cc \
    voro/wall.cc \

HEADERS += \
    src/box.h \
    src/cell.h \
    src/config.h \
    src/contact.h \
    src/event.h \
    src/grain.h \
    src/grid.h \
    src/in_out.h \
    src/main.h \
    src/mesh.h \
    src/node.h \
    src/parameter.h \
    src/profile.h \
    src/run.h \
    src/run_post_process.h \
    src/run_visu.h \
    src/slice.h \
    src/vector.h \
    voro/c_loops.hh \
    voro/cell.hh \
    voro/common.hh \
    voro/config.hh \
    voro/container.hh \
    voro/container_prd.hh \
    voro/pre_container.hh \
    voro/rad_option.hh \
    voro/unitcell.hh \
    voro/v_base.hh \
    voro/v_compute.hh \
    voro/voro++.hh \
    voro/wall.hh \
    voro/worklist.hh
