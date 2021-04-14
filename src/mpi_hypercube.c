/*
 * mpi_hypercube -- Implements a hypercube-interconnect network topology using
 * OpenMPi Copyright (c) 2021 Ángel Pérez <angel@ttm.sh>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define PROGNAME "mpi_hypercube"
#define DISTRIB_RANK 0
#define TAG_FINAL_RESULT 42

#define max(a, b)                                                             \
    __extension__({                                                           \
        __typeof(a) _a = a;                                                   \
        __typeof(b) _b = b;                                                   \
        _a > _b ? _a : _b;                                                    \
    })
#define MPI_Check(v)                                                          \
    __extension__({                                                           \
        __typeof(v) _v = v;                                                   \
        if (_v != MPI_SUCCESS)                                                \
            handle_error(_v, #v);                                             \
    })
#define logf(f, ...) printf(PROGNAME "(%d): " f "\n", g_rank, ##__VA_ARGS__)

static int g_rank = -1, g_size = -1;

/* Generic MPI error handler.
 * This function gets called from within the MPI_Check() macro in case a MPI
 * call does not succeed. Do not attempt to call this handler manually.
 * This function does not return.
 * @mpi_error: MPI status code
 * @expr: Failing expression
 */
static void __attribute__((noreturn))
handle_error(int mpi_error, const char *expr)
{
    char msg_buf[BUFSIZ];
    int msg_len = -1;

    if (MPI_Error_string(mpi_error, msg_buf, &msg_len) != MPI_SUCCESS
        && msg_len > 0) {
        fprintf(stderr, PROGNAME "(%d): MPI error %d (`%s')\n", g_rank,
            mpi_error, expr);
    } else {
        fprintf(stderr, PROGNAME "(%d): MPI error %d (`%s'): %s\n", g_rank,
            mpi_error, expr, msg_buf);
    }

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    MPI_Finalize();
    _exit(EXIT_FAILURE);
}

/* Parse the dimensions of the hypercube from the command line arguments.
 * Returns -1 on failure, parsed value on success
 * @str: String to be parsed
 */
static int parse_dimensions(char *str)
{
    errno = 0;
    char *endptr;
    long result = strtol(str, &endptr, 10);

    if (endptr == str)
        return -1; /* nothing parsed */

    if ((result == INT_MAX || result == INT_MIN) && errno == ERANGE)
        return -1; /* out of range */
    return (int)result;
}

/* Returns 1 if the given character is part of a valid floating-point numeric
 * entity.
 * @c: Character to be checked
 */
static inline int __attribute__((always_inline)) is_numeric_char(char c)
{
    return isdigit(c) || c == '.' || c == '-';
}

/* Reads the input file and scans for numeric entities and sends out these
 * values across the network for every peer to process.
 * @path: Path to the file containing the data
 */
static void perform_distribution(const char *path)
{
    /* Open input file. */
    errno = 0;
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr,
            PROGNAME "(%d): error: could not open file `%s' for reading: %s\n",
            g_rank, path, strerror(errno));
        MPI_Check(MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE));
        MPI_Check(MPI_Finalize());
        _exit(EXIT_FAILURE);
    }

    /* Keep reading from file until we're done. */
    int n = 0;
    char c = 0, item_buf[BUFSIZ];
    for (;;) {
        size_t i = 0;

        /* Copy characters into our item buffer as long as they are part of a
         * numeric entity. */
        while (i < BUFSIZ && is_numeric_char((c = (char)fgetc(fp))))
            item_buf[i++] = c;
        item_buf[i] = '\0';

        if (c == EOF) /* EOF reached -- we're good to go */
            break;

        if (n >= g_size - 1) {
            fprintf(stderr,
                PROGNAME "(%d): warning: too many numeric "
                         "entities on the list. %d values were "
                         "expected but this is the %dth. Only first %d "
                         "values will be taken into account\n",
                g_rank, g_size - 1, 1 + n, g_size - 1);
            break;
        }

        if (i >= BUFSIZ - 1) {
            /* Buffer overflow -- skip entity */
            logf("warning: skipping entity overflowing buffer");
            continue;
        }

        /* Parse double (see <https://stackoverflow.com/a/5581058>) */
        char *endptr;
        double val = strtod(item_buf, &endptr);
        if (*endptr != '\0') {
            logf("warning: skipping invalid entity (`%s')", item_buf);
            continue;
        }

        /* Send out this value to the rest of processes
         * we send it to (1 + n) because first worker is always the
         * distributor process. */
        MPI_Check(MPI_Bsend(&val, 1, MPI_DOUBLE, 1 + n, 0, MPI_COMM_WORLD));
        n++;
    }

    if (n < g_size - 1) {
        fprintf(stderr,
            PROGNAME "(%d): error: invalid number of values on the list. "
                     "Expected exactly %d but got %d.\n",
            g_rank, g_size - 1, n);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        MPI_Finalize();
        _exit(EXIT_FAILURE);
    }

    fclose(fp);

    /* Receive the minimum value back from the worker processes. */
    MPI_Status status;
    double minimum_value;
    MPI_Recv(&minimum_value, 1, MPI_DOUBLE, MPI_ANY_SOURCE, TAG_FINAL_RESULT,
        MPI_COMM_WORLD, &status);
    MPI_Check(status.MPI_ERROR);
    printf("%lf\n", minimum_value);
}

static void get_neighbors(int *neighbors, int dim)
{
    for (int i = 0; i < dim; i++)
        neighbors[i] = 1 + ((g_rank - 1) ^ (1 << i));
}

static void do_work(int dim)
{
    MPI_Status status;

    /* Receive value from the distributor process. */
    double distrib_val;
    MPI_Recv(&distrib_val, 1, MPI_DOUBLE, DISTRIB_RANK, MPI_ANY_TAG,
        MPI_COMM_WORLD, &status);
    MPI_Check(status.MPI_ERROR);

    /* Obtain neighbor processes. */
    int neighbors[dim];
    get_neighbors(neighbors, dim);

    /* Receive value from neighbor and calculate maximum. */
    double neighbor_val;
    for (int i = 0; i < dim; i++) {
        MPI_Check(MPI_Bsend(
            &distrib_val, 1, MPI_DOUBLE, neighbors[i], 0, MPI_COMM_WORLD));
        MPI_Recv(&neighbor_val, 1, MPI_DOUBLE, neighbors[i], MPI_ANY_TAG,
            MPI_COMM_WORLD, &status);
        MPI_Check(status.MPI_ERROR);
        distrib_val = max(distrib_val, neighbor_val);
    }

    /* Send out the minimum value back to the distributor process. */
    MPI_Check(MPI_Bsend(&distrib_val, 1, MPI_DOUBLE, DISTRIB_RANK,
        TAG_FINAL_RESULT, MPI_COMM_WORLD));
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        printf("usage: mpi_hypercube DIMENSION INPUT_FILE\n\n");
        return EXIT_SUCCESS;
    }

    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fprintf(stderr, PROGNAME ": error: MPI initialization failed\n");
        return EXIT_FAILURE;
    }

    MPI_Check(MPI_Comm_rank(MPI_COMM_WORLD, &g_rank));
    MPI_Check(MPI_Comm_size(MPI_COMM_WORLD, &g_size));

    /* Parse and check dimension for the hypercube topology. */
    int dim = parse_dimensions(argv[1]);
    if (dim < 2) {
        fprintf(stderr, PROGNAME "(%d): error: invalid dimension (%d)\n",
            g_rank, dim);
        MPI_Check(MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE));
        MPI_Check(MPI_Finalize());
        return EXIT_FAILURE;
    }

    /* Are there enough processes for this topology? */
    int num_expected_slots = 1 + (int)powl(2L, dim);
    if (g_size < num_expected_slots) {
        fprintf(stderr,
            PROGNAME "(%d): error: no enough slots for toroid topology. Got "
                     "%d when %d processes were expected.\n",
            g_rank, g_size, num_expected_slots);
        MPI_Check(MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE));
        MPI_Check(MPI_Finalize());
        return EXIT_FAILURE;
    }

    /* Is this the distributor process? */
    if (g_rank == DISTRIB_RANK)
        perform_distribution(argv[2]);
    else
        do_work(dim);

    MPI_Check(MPI_Finalize());
    return EXIT_SUCCESS;
}
