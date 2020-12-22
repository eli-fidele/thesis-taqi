# Swing Stacks: Predicting political party affiliation via model stacks

This repository contains the final project for Math 243: Statisical Learning at Reed College in Fall 2020. Authors of this project are Shisham Adhikari, Maggie Slein, and Grayson White. The final paper can be found in the `technical-report.pdf` file and the final presentation slides in the `final-presentation.pdf` file. Below is the abstract.

## Abstract

As American politics has become increasingly polarized over the last several decades, predicting the likelihood of presidential victories has become more difficult. With just a few states' electoral votes deciding the figurehead of the executive branch, two out of the four last elections have been decided based on the electoral college and not on the popular vote. Given this, many previous models used to predict election outcomes have failed to predict the true winner. We aim to address how socio-political factors influence political party affiliation in two ways: First, by using the `survey` package to create the "classic" logistic regression models with weights, and then by using `tidymodels` and `stacks` to layer multiple models and model types to better predict party affiliation using several predictors from the General Social Survey (GSS) dataset from 2000-2016. 
