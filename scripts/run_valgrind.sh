# This script runs a full Valgrind memory check of the application given as the only argument and saves the output to 'valgrind.out'.
#!/bin/bash
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind.out $1
