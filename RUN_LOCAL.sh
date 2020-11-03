#!/bin/bash

for value in {104..875}
do
	sage structure_distinguished.sage $value
done

echo DONE!
