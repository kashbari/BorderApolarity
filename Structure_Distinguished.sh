#!/bin/bash

for value in {104..526}
do
	sage structure_distinguished.sage $value
done

echo DONE!
