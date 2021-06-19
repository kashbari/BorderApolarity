#!/bin/bash

for value in {0..607}
do
	sage structure3_16_flag.sage $value
done

echo DONE!
