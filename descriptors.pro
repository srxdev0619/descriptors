QT += core
QT -= gui

CONFIG += c++11

TARGET = descriptors
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    descriptor3.cpp \
    checker.cpp \
    common_structs.cpp \
    symmetric_functions.cpp

HEADERS += \
    descriptor3.h \
    checker.h \
    common_structs.h \
    symmetric_functions.h
