#!/bin/bash

for value in {0..526}
do
	sage structure_distinguished.sage $value
done

echo DONE!
