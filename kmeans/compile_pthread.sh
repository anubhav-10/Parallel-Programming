#!/bin/bash
g++ -Ofast -o pthread main_pthread.c lab1_io.c pthread_v2.cpp -lpthread -lgomp

