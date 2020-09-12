#!/bin/bash

for value in {36..875}
do
	sage structure_distinguished.sage $value
done

echo DONE!
