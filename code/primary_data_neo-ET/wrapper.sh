#!/bin/bash

#copy all required data locally
#data.preparation.sh

#Run Sleuth and initial enrichments
R --vanilla < sleuth.run.R
R --vanilla < sleuth.gsea.R

#Run downstream analyses (new neo-adjuvant data)
R --vanilla < downstream.neo.treatment.basic.R
R --vanilla < downstream.neo.treatment.regulators.R
R --vanilla < downstream.neo.response.all.R

#Run downstream analyses (pervasiveness)
R --vanilla < downstream.pervasiveness.basic.R
R --vanilla < downstream.pervasiveness.fgsea.R
R --vanilla < downstream.pervasiveness.string.R
