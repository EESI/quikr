#! /bin/sh -v
mllex substitution.lex
mlyacc substitution.grm
mlton -link-opt -lJudy -link-opt -lz score.mlb gzip.c
