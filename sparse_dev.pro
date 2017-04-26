TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -pthread

SOURCES += main.cpp \
    src/lib/sparse_matrix/msr_matrix.cpp \
    src/lib/sparse_matrix/msr_algo.cpp \
    src/workers/solver.cpp \
    src/lib/threads/thread_vector_utils.cpp \
    src/lib/threads/thread_handler.cpp \
    src/lib/sparse_matrix/msr_thread_handler.cpp \
    src/lib/sparse_matrix/msr_thread_dqgmres_solver.cpp

HEADERS += \
    src/lib/sparse_matrix/msr_matrix.h \
    src/lib/sparse_matrix/msr_algo.h \
    src/workers/solver.h \
    src/lib/containers/cycle_buf.h \
    src/lib/threads/thread_vector_utils.h \
    src/lib/threads/thread_handler.h \
    src/lib/sparse_matrix/msr_thread_handler.h \
    src/lib/sparse_matrix/msr_thread_dqgmres_solver.h

INCLUDEPATH += src/
INCLUDEPATH += src/lib

