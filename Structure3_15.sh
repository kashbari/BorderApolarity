#!/bin/bash

for value in {0..527}
do
	sage structure3_15.sage $value
done

echo DONE!
