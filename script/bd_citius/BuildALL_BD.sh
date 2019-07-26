#!/usr/bin/env bash

rm -Rf /media/glusterfs/genomics/Databases/AFS20_4/*
rm -Rf /media/glusterfs/genomics/Databases/AFS20_8/*
rm -Rf /media/glusterfs/genomics/Databases/AFS20_16/*
rm -Rf /media/glusterfs/genomics/Databases/AFS20_32/*
rm -Rf /media/glusterfs/genomics/Databases/AFS20_64/*
rm -Rf /media/glusterfs/genomics/Databases/AFS20_128/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31_4/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31_8/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31_16/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31_32/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31_64/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31_128/*
rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_4/*
rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_8/*
rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_16/*
rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_32/*
rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_64/*
rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_128/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_4/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_8/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_16/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_32/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_64/*
rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_128/*


for i in 1 2
do
    ./BuildGeneric_BD.sh AFS20 4
    ./BuildGeneric_BD.sh AFS20 8
    ./BuildGeneric_BD.sh AFS20 16
    ./BuildGeneric_BD.sh AFS20 32
    ./BuildGeneric_BD.sh AFS20 64
    ./BuildGeneric_BD.sh AFS20 128
    ./BuildGeneric_BD.sh AFS31 4
    ./BuildGeneric_BD.sh AFS31 8
    ./BuildGeneric_BD.sh AFS31 16
    ./BuildGeneric_BD.sh AFS31 32
    ./BuildGeneric_BD.sh AFS31 64
    ./BuildGeneric_BD.sh AFS31 128
    ./BuildGeneric_BD.sh AFS20RS90 4
    ./BuildGeneric_BD.sh AFS20RS90 8
    ./BuildGeneric_BD.sh AFS20RS90 16
    ./BuildGeneric_BD.sh AFS20RS90 32
    ./BuildGeneric_BD.sh AFS20RS90 64
    ./BuildGeneric_BD.sh AFS20RS90 128
    ./BuildGeneric_BD.sh AFS31RS90 4
    ./BuildGeneric_BD.sh AFS31RS90 8
    ./BuildGeneric_BD.sh AFS31RS90 16
    ./BuildGeneric_BD.sh AFS31RS90 32
    ./BuildGeneric_BD.sh AFS31RS90 64
    ./BuildGeneric_BD.sh AFS31RS90 128

    rm -Rf /media/glusterfs/genomics/Databases/AFS20_4/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS20_8/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS20_16/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS20_32/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS20_64/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS20_128/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31_4/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31_8/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31_16/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31_32/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31_64/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31_128/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_4/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_8/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_16/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_32/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_64/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS20RS90_128/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_4/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_8/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_16/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_32/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_64/*
    rm -Rf /media/glusterfs/genomics/Databases/AFS31RS90_128/*

done

./BuildGeneric_BD.sh AFS20 4
./BuildGeneric_BD.sh AFS20 8
./BuildGeneric_BD.sh AFS20 16
./BuildGeneric_BD.sh AFS20 32
./BuildGeneric_BD.sh AFS20 64
./BuildGeneric_BD.sh AFS20 128
./BuildGeneric_BD.sh AFS31 4
./BuildGeneric_BD.sh AFS31 8
./BuildGeneric_BD.sh AFS31 16
./BuildGeneric_BD.sh AFS31 32
./BuildGeneric_BD.sh AFS31 64
./BuildGeneric_BD.sh AFS31 128
./BuildGeneric_BD.sh AFS20RS90 4
./BuildGeneric_BD.sh AFS20RS90 8
./BuildGeneric_BD.sh AFS20RS90 16
./BuildGeneric_BD.sh AFS20RS90 32
./BuildGeneric_BD.sh AFS20RS90 64
./BuildGeneric_BD.sh AFS20RS90 128
./BuildGeneric_BD.sh AFS31RS90 4
./BuildGeneric_BD.sh AFS31RS90 8
./BuildGeneric_BD.sh AFS31RS90 16
./BuildGeneric_BD.sh AFS31RS90 32
./BuildGeneric_BD.sh AFS31RS90 64
./BuildGeneric_BD.sh AFS31RS90 128