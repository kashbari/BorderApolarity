#!/bin/bash

for value in {522..526}
do
	sage structure3_15.sage $value
done

echo DONE!
