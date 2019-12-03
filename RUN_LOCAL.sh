#!/bin/bash

for value in {302..305}
do
	sage structure.sage $value
done

echo DONE!
