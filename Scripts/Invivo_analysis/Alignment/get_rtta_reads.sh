#!/bin/bash

module load samtools

samtools view -F 4 $1 | wc -l > $1.rtta_reads
