#!/bin/bash

if [[ $1 == "" ]]
then
	echo " * -------------------------------- *"
	echo " * OptiMass model src & lib builder *"
	echo " * -------------------------------- *"
	echo " Usage: $./build_model.sh <your_model.xml>"
	echo " "
	exit -1;
fi

modelfile=$1
modelcard="./model/model_card.xml"

if [ -f $modelfile ]
then
    echo ""
    echo " *processing ${modelfile}..."
    cp $modelfile $modelcard;
    ./build_model_dictionary;
    make model;
else
    "There's no modelfile named ${modelfile}"
fi

exit 0;

