#!/bin/bash

PERF_DATA_DIR="data/perf_data"
mkdir -p $PERF_DATA_DIR

echo "Running performance profile for bc_model"
Rscript -e "
  library(dde)
  Rprof('$PERF_DATA_DIR/bc_perf.out', interval = 0.001)
  source('Breast_Cancer_Model/bc_model.R')
  Rprof(NULL)
  prof <- summaryRprof('$PERF_DATA_DIR/bc_perf.out')
  write.csv(prof\$by.self, '$PERF_DATA_DIR/bc_perf_report.csv')
  cat('Profile saved to $PERF_DATA_DIR/bc_perf_report.csv\n')
"
if [ $? -eq 0 ]; then
  echo "Breast Cancer performance test passed"
else
  echo "Breast Cancer performance test failed"
  exit 1
fi

echo "Running performance profile for cv_model"
Rscript -e "
  library(dde)
  Rprof('$PERF_DATA_DIR/cv_perf.out', interval = 0.001)
  source('Cardiovascular_Model/cv_model.R')
  Rprof(NULL)
  prof <- summaryRprof('$PERF_DATA_DIR/cv_perf.out')
  write.csv(prof\$by.self, '$PERF_DATA_DIR/cv_perf_report.csv')
  cat('Profile saved to $PERF_DATA_DIR/cv_perf_report.csv\n')
"
if [ $? -eq 0 ]; then
  echo "Cardiovascular performance test passed"
else
  echo "Cardiovascular performance test failed"
  exit 1
fi

echo "Running performance profile for cw_model"
Rscript -e "
  library(dde)
  Rprof('$PERF_DATA_DIR/cw_perf.out', interval = 0.001)
  source('Cobweb_Model/cw_model.R')
  Rprof(NULL)
  prof <- summaryRprof('$PERF_DATA_DIR/cw_perf.out')
  write.csv(prof\$by.self, '$PERF_DATA_DIR/cw_perf_report.csv')
  cat('Profile saved to $PERF_DATA_DIR/cw_perf_report.csv\n')
"
if [ $? -eq 0 ]; then
  echo "Cobweb performance test passed"
else
  echo "Cobweb performance test failed"
  exit 1
fi

echo "Running performance profile for gi_model"
Rscript -e "
  library(dde)
  Rprof('$PERF_DATA_DIR/gi_perf.out', interval = 0.001)
  source('Glucose_Insulin_Model/gi_model.R')
  Rprof(NULL)
  prof <- summaryRprof('$PERF_DATA_DIR/gi_perf.out')
  write.csv(prof\$by.self, '$PERF_DATA_DIR/gi_perf_report.csv')
  cat('Profile saved to $PERF_DATA_DIR/gi_perf_report.csv\n')
"
if [ $? -eq 0 ]; then
  echo "Glucose Insulin performance test passed"
else
  echo "Glucose Insulin performance test failed"
  exit 1
fi