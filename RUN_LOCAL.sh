#!/bin/bash

for value in {253..305}
do
	sage structure.sage $value
done

echo DONE!
