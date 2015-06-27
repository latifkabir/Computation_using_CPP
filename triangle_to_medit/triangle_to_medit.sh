#!/bin/bash
#
g++ -c triangle_to_medit.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_to_medit.cpp"
  exit
fi
#
g++ triangle_to_medit.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_to_medit.o."
  exit
fi
rm triangle_to_medit.o
#
mv a.out ~/bincpp/$ARCH/triangle_to_medit
#
echo "Executable installed as ~/bincpp/$ARCH/triangle_to_medit"
