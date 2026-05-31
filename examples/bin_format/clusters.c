#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct Grid {
    float dx, dy, dz;
    float ox, oy, oz;
    uint32_t nx, ny, nz;
    uint64_t npoints;
    float* data;
} Grid;

// -----------------------------------------------------------------------------
static int system_is_little_endian(void) {
    uint16_t x = 1;
    return (*(uint8_t*)&x) == 1;
}

// -----------------------------------------------------------------------------
static void byteswap_u32_array(uint32_t* buf, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        uint32_t v = buf[i];
        buf[i] =
            ((v & 0x000000FFu) << 24) |
            ((v & 0x0000FF00u) << 8)  |
            ((v & 0x00FF0000u) >> 8)  |
            ((v & 0xFF000000u) >> 24);
    }
}

// -----------------------------------------------------------------------------
static uint32_t read_u32_le(FILE* file) {
    unsigned char b[4];
    if (fread(b, 1, 4, file) != 4) {
        fprintf(stderr, "Unexpected EOF reading uint32\n");
        exit(EXIT_FAILURE);
    }
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) | ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

// -----------------------------------------------------------------------------
static float read_f32_le(FILE* file) {
    uint32_t u = read_u32_le(file);
    float v;
    memcpy(&v, &u, sizeof(float));
    return v;
}

// -----------------------------------------------------------------------------
static void write_u32_le(FILE* file, uint32_t v) {
    unsigned char b[4];
    b[0] =  v        & 0xFF;
    b[1] = (v >> 8 ) & 0xFF;
    b[2] = (v >> 16) & 0xFF;
    b[3] = (v >> 24) & 0xFF;
    fwrite(b, 1, 4, file);
}

// -----------------------------------------------------------------------------
static void write_f32_le(FILE* file, float v) {
uint32_t u;
    memcpy(&u, &v, sizeof(float));
    write_u32_le(file, u);
}

// -----------------------------------------------------------------------------
static int zeros_like(Grid* new_grid, Grid* ref) {
    new_grid->nx = ref->nx;
    new_grid->ny = ref->ny;
    new_grid->nz = ref->nz;
    new_grid->dx = ref->dx;
    new_grid->dy = ref->dy;
    new_grid->dz = ref->dz;
    new_grid->ox = ref->ox;
    new_grid->oy = ref->oy;
    new_grid->oz = ref->oz;
    new_grid->npoints = ref->npoints;

    new_grid->data = (float*)malloc(ref->npoints * sizeof(float));
    if (!new_grid->data) { fprintf(stderr, "malloc failed\n"); return 1; }
    for (uint64_t i = 0; i < ref->npoints; ++i) new_grid->data[i] = 0;

    return 0;
}

// -----------------------------------------------------------------------------
static int read_bin(const char* path_bin, Grid* grid) {
    FILE* file = fopen(path_bin, "rb");
    if (!file) { perror("fopen input"); return 1; }

    grid->nx = read_u32_le(file);
    grid->ny = read_u32_le(file);
    grid->nz = read_u32_le(file);
    grid->dx = read_f32_le(file);
    grid->dy = read_f32_le(file);
    grid->dz = read_f32_le(file);
    grid->ox = read_f32_le(file);
    grid->oy = read_f32_le(file);
    grid->oz = read_f32_le(file);
    grid->npoints = (uint64_t)grid->nx * (uint64_t)grid->ny * (uint64_t)grid->nz;

    grid->data = (float*)malloc(grid->npoints * sizeof(float));
    if (!grid->data) { fprintf(stderr, "malloc failed\n"); fclose(file); return 1; }

    if (system_is_little_endian()) {
        size_t got = fread(grid->data, sizeof(float), (size_t)grid->npoints, file);
        if (got != (size_t)grid->npoints) { fprintf(stderr, "Unexpected EOF reading data\n"); free(grid->data); fclose(file); return 1; }
        fclose(file);
        return 0;
    }

    //// read as u32 and byteswap to host order, then memcpy into float array
    uint32_t* tmp = (uint32_t*)malloc(grid->npoints * sizeof(uint32_t));
    if (!tmp) { fprintf(stderr, "malloc tmp failed\n"); free(grid->data); fclose(file); return 1; }
    size_t got = fread(tmp, sizeof(uint32_t), (size_t)grid->npoints, file);
    if (got != (size_t)grid->npoints) { fprintf(stderr, "Unexpected EOF reading data\n"); free(tmp); free(grid->data); fclose(file); return 1; }
    byteswap_u32_array(tmp, (size_t)grid->npoints);
    memcpy(grid->data, tmp, grid->npoints * sizeof(float));
    free(tmp);

    fclose(file);
    return 0;
}

// -----------------------------------------------------------------------------
static int write_bin(const char* path_bin, Grid* grid) {
    FILE* file = fopen(path_bin, "wb");
    if (!file) { perror("fopen output"); free(grid->data); return 1; }

    write_u32_le(file, grid->nx);
    write_u32_le(file, grid->ny);
    write_u32_le(file, grid->nz);
    write_f32_le(file, grid->dx);
    write_f32_le(file, grid->dy);
    write_f32_le(file, grid->dz);
    write_f32_le(file, grid->ox);
    write_f32_le(file, grid->oy);
    write_f32_le(file, grid->oz);

    int little2 = system_is_little_endian();
    if (little2) {
        size_t wrote = fwrite(grid->data, sizeof(float), (size_t)grid->npoints, file);
        if (wrote != (size_t)grid->npoints) { fprintf(stderr, "Error writing data\n"); free(grid->data); fclose(file); return 1; }
        fclose(file);
        return 0;
    }

    uint32_t* tmp = (uint32_t*)malloc(grid->npoints * sizeof(uint32_t));
    if (!tmp) { fprintf(stderr, "malloc tmp failed\n"); free(grid->data); fclose(file); return 1; }
    memcpy(tmp, grid->data, grid->npoints * sizeof(float));
    byteswap_u32_array(tmp, (size_t)grid->npoints);
    size_t wrote = fwrite(tmp, sizeof(uint32_t), (size_t)grid->npoints, file);
    free(tmp);
    if (wrote != (size_t)grid->npoints) { fprintf(stderr, "Error writing data\n"); free(grid->data); fclose(file); return 1; }

    fclose(file);
    return 0;
}

// -----------------------------------------------------------------------------
int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s input.bin output.bin threshold\n", argv[0]);
        return 1;
    }
    const char* path_in  = argv[1];
    const char* path_out = argv[2];
    float thr = (float)atof(argv[3]);
    int status;

    Grid grid = { 0 };
    status = read_bin(path_in, &grid);
    if (status != 0) { return status; }

    Grid clusters = { 0 };
    status = zeros_like(&clusters, &grid);
    if (status != 0) { return status; }


    uint32_t nz = grid.nz;
    uint32_t nyz = grid.ny * nz;
    float current_cluster = 1.0f;
    // [WIP] sample logic, this won't actually provide a proper segmentation
    for (uint64_t i = 0; i < grid.npoints; ++i) {
        float point = grid.data[i];
        int idx_xn = i - nyz;
        int idx_yn = i % nyz >= nz ? i - nz : -1;
        int idx_zn = i % nz > 0.0f ? i - 1 : -1;

        if (grid.data[i] < thr) continue;
        if (idx_xn >= 0 && clusters.data[idx_xn]) { clusters.data[i] = clusters.data[idx_xn]; continue; }
        if (idx_yn >= 0 && clusters.data[idx_yn]) { clusters.data[i] = clusters.data[idx_yn]; continue; }
        if (idx_zn >= 0 && clusters.data[idx_zn]) { clusters.data[i] = clusters.data[idx_zn]; continue; }

        clusters.data[i] = current_cluster++;
    }

    status = write_bin(path_out, &clusters);
    if (status != 0) { return status; }

    free(grid.data);
    free(clusters.data);

    printf("Wrote %s (nx=%u ny=%u nz=%u) thresholded at %g\n", path_out, grid.nx, grid.ny, grid.nz, (double)thr);
    return 0;
}

// -----------------------------------------------------------------------------
