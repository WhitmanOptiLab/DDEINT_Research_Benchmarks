#!/bin/bash

echo "Running benchmark for bc_model"
Rscript Breast_Cancer_Model/bc_model.R
if [ $? -eq 0 ]; then
  echo "Breast Cancer benchmark passed"
else
  echo "Breast Cancer benchmark failed"
  exit 1
fi

echo "Running benchmark for cv_model"
Rscript Cardiovascular_Model/cv_model.R
if [ $? -eq 0 ]; then
  echo "Cardiovascular benchmark passed"
else
  echo "Cardiovascular benchmark failed"
  exit 1
fi

echo "Running benchmark for cw_model"
Rscript Cobweb_Model/cw_model.R
if [ $? -eq 0 ]; then
  echo "Cobweb benchmark passed"
else
  echo "Cobweb benchmark failed"
  exit 1
fi

echo "Running benchmark for gi_model"
Rscript Glucose_Insulin_Model/gi_model.R
if [ $? -eq 0 ]; then
  echo "Glucose Insulin benchmark passed"
else
  echo "Glucose Insulin benchmark failed"
  exit 1
fi

echo "All dde R benchmarks complete. Results saved to data"