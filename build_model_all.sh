#!/bin/bash

for file_xml in ./model/example_models/*.xml
do
    echo $file_xml
    cp $file_xml ./model/model_card.xml
    ./build_model_dictionary

done

exit 0

