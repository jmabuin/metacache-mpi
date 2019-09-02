#!/usr/bin/env bash

#rm -Rf /home/josemanuel.abuin/Genomics/Databases/AFS31RS90_8/*
#rm -Rf /home/josemanuel.abuin/Genomics/Databases/AFS31RS90_16/*
#rm -Rf /home/josemanuel.abuin/Genomics/Databases/AFS31RS90_32/*
rm -Rf /home/josemanuel.abuin/Genomics/Databases/AFS31RS90_64/*
rm -Rf /home/josemanuel.abuin/Genomics/Databases/AFS31RS90_128/*


for i in 1 2
do
 #   ./BuildGeneric_Local_BD.sh AFS31RS90 8
 #   ./BuildGeneric_Local_BD.sh AFS31RS90 16
 #   ./BuildGeneric_Local_BD.sh AFS31RS90 32
    ./BuildGeneric_Local_BD.sh AFS31RS90 64
    ./BuildGeneric_Local_BD.sh AFS31RS90 128

  #  rm -Rf /home/josemanuel.abuin/Genomics/Databases/AFS31RS90_8/*
  #  rm -Rf /home/josemanuel.abuin/Genomics/Databases/AFS31RS90_16/*
  #  rm -Rf /home/josemanuel.abuin/Genomics/Databases/AFS31RS90_32/*
    rm -Rf /home/josemanuel.abuin/Genomics/Databases/AFS31RS90_64/*
    rm -Rf /home/josemanuel.abuin/Genomics/Databases/AFS31RS90_128/*

done

#./BuildGeneric_Local_BD.sh AFS31RS90 8
#./BuildGeneric_Local_BD.sh AFS31RS90 16
#./BuildGeneric_Local_BD.sh AFS31RS90 32
./BuildGeneric_Local_BD.sh AFS31RS90 64
./BuildGeneric_Local_BD.sh AFS31RS90 128