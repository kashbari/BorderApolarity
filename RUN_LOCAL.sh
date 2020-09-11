#!/bin/bash

for value in {14..875}
do
	sage structure_distinguished.sage $value
done

echo DONE!
