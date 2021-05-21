#!/bin/bash


easythresh_conj.sh zstat2 zstat3 mask 2.6 0.05 mean_func conj23_pos
easythresh_conj.sh thresh_conj23_pos zstat4 mask 2.6 0.05 mean_func conj234_pos
easythresh_conj.sh thresh_conj234_pos zstat5 mask 2.6 0.05 mean_func conj2345_pos

easythresh_conj.sh zstat2_neg zstat3_neg mask 2.6 0.05 mean_func conj23_neg
easythresh_conj.sh thresh_conj23_neg zstat4_neg mask 2.6 0.05 mean_func conj234_neg
easythresh_conj.sh thresh_conj234_neg zstat5_neg mask 2.6 0.05 mean_func conj2345_neg

fslmaths thresh_conj2345_neg -mul -1 -add thresh_conj2345_pos thresh_conj2345

