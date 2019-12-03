#!/bin/bash

for value in {131..527}
do
	sage structure3_15.sage $value
done

echo DONE!
