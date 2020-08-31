#!/bin/bash

for value in {2,3,4}
do
	sage structure4_16.sage $value
done

echo DONE!
