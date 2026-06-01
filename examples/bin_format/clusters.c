#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

typedef struct Grid {
    float dx, dy, dz;
    float ox, oy, oz;
    uint32_t nx, ny, nz;
    uint64_t npoints;
    float* data;
} Grid;

typedef struct Node {
    uint64_t value;
    struct Node* next;
} Node;

// -----------------------------------------------------------------------------
static bool system_is_little_endian(void) {
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

    if (system_is_little_endian()) {
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
static Node* extend_linked_list(
    Node* graph, // start of the linked list (i.e. "graph" array)
    Node* node,  // pointer to the current node's index (i.e. "graph[i]")
    int idx_neighbor // flattened index for the arrays, -1 if out of bounds
) {
    if (idx_neighbor < 0) return node;

    Node* neighbor = graph + idx_neighbor;

    if (neighbor->next) return node; // already in the list, skip
    neighbor->value = (uint64_t)idx_neighbor;

    neighbor->next = node->next;
    node->next = neighbor;

    return node->next;
}

// -----------------------------------------------------------------------------
static Node* pop_linked_list(Node* node) {
    if (!node) return NULL;
    Node* next = node->next;
    node->next = NULL; // mark as visited
    return next;
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

    Grid smif = { 0 }; // input grid
    status = read_bin(path_in, &smif);
    if (status != 0) { return status; }

    Grid clusters = { 0 }; // output grid
    status = zeros_like(&clusters, &smif);
    if (status != 0) { return status; }

    // linked list for traversing grid points' neighbors
    Node* graph = (Node*)malloc(smif.npoints * sizeof(Node));
    if (!graph) { fprintf(stderr, "malloc failed\n"); return 1; }
    for (uint64_t i = 0; i < smif.npoints; ++i) {
        graph[i].value = 0; graph[i].next = NULL;
    }

    uint32_t nx  = smif.nx;
    uint32_t ny  = smif.ny;
    uint32_t nz  = smif.nz;
    uint32_t nyz = smif.ny * nz;
    uint64_t npoints = smif.npoints;
    float current_cluster = 0.0f; // float (despite cluster labels being integers) because the output grid is float

    for (uint64_t i = 0; i < npoints; ++i) {
        if (smif.data[i] < thr || clusters.data[i]) continue;

        current_cluster++;

        Node* head = graph + i;
        head->value = (uint64_t)i;

        while (head) {
            uint64_t idx = head->value;

            if (smif.data[idx] < thr || clusters.data[idx]) {
                head = pop_linked_list(head);
                continue;
            }

            clusters.data[idx] = current_cluster;

            int x =  idx / nyz;
            int y = (idx % nyz) / nz;
            int z =  idx % nz;

            bool xmin = x == 0;
            bool xmax = x == (int)nx - 1;
            bool ymin = y == 0;
            bool ymax = y == (int)ny - 1;
            bool zmin = z == 0;
            bool zmax = z == (int)nz - 1;

            int x0y_z_ = xmin?-1: (x-1) * nyz +  y    * nz +  z   ;
            int x1y_z_ = xmax?-1: (x+1) * nyz +  y    * nz +  z   ;
            int x_y0z_ = ymin?-1:  x    * nyz + (y-1) * nz +  z   ;
            int x_y1z_ = ymax?-1:  x    * nyz + (y+1) * nz +  z   ;
            int x_y_z0 = zmin?-1:  x    * nyz +  y    * nz + (z-1);
            int x_y_z1 = zmax?-1:  x    * nyz +  y    * nz + (z+1);

            int x0y0z_ = (xmin || ymin)?-1: (x-1) * nyz + (y-1) * nz +  z   ;
            int x0y1z_ = (xmin || ymax)?-1: (x-1) * nyz + (y+1) * nz +  z   ;
            int x0y_z0 = (xmin || zmin)?-1: (x-1) * nyz +  y    * nz + (z-1);
            int x0y_z1 = (xmin || zmax)?-1: (x-1) * nyz +  y    * nz + (z+1);

            int x1y0z_ = (xmax || ymin)?-1: (x+1) * nyz + (y-1) * nz +  z   ;
            int x1y1z_ = (xmax || ymax)?-1: (x+1) * nyz + (y+1) * nz +  z   ;
            int x1y_z0 = (xmax || zmin)?-1: (x+1) * nyz +  y    * nz + (z-1);
            int x1y_z1 = (xmax || zmax)?-1: (x+1) * nyz +  y    * nz + (z+1);

            int x_y0z0 = (ymin || zmin)?-1:  x    * nyz + (y-1) * nz + (z-1);
            int x_y0z1 = (ymin || zmax)?-1:  x    * nyz + (y-1) * nz + (z+1);
            int x_y1z0 = (ymax || zmin)?-1:  x    * nyz + (y+1) * nz + (z-1);
            int x_y1z1 = (ymax || zmax)?-1:  x    * nyz + (y+1) * nz + (z+1);

            int x0y0z0 = (xmin || ymin || zmin)?-1: (x-1) * nyz + (y-1) * nz + (z-1);
            int x0y0z1 = (xmin || ymin || zmax)?-1: (x-1) * nyz + (y-1) * nz + (z+1);
            int x0y1z0 = (xmin || ymax || zmin)?-1: (x-1) * nyz + (y+1) * nz + (z-1);
            int x0y1z1 = (xmin || ymax || zmax)?-1: (x-1) * nyz + (y+1) * nz + (z+1);

            int x1y0z0 = (xmax || ymin || zmin)?-1: (x+1) * nyz + (y-1) * nz + (z-1);
            int x1y0z1 = (xmax || ymin || zmax)?-1: (x+1) * nyz + (y-1) * nz + (z+1);
            int x1y1z0 = (xmax || ymax || zmin)?-1: (x+1) * nyz + (y+1) * nz + (z-1);
            int x1y1z1 = (xmax || ymax || zmax)?-1: (x+1) * nyz + (y+1) * nz + (z+1);

            extend_linked_list(graph, head, x0y_z_);
            extend_linked_list(graph, head, x1y_z_);
            extend_linked_list(graph, head, x_y0z_);
            extend_linked_list(graph, head, x_y1z_);
            extend_linked_list(graph, head, x_y_z0);
            extend_linked_list(graph, head, x_y_z1);
            extend_linked_list(graph, head, x0y0z_);
            extend_linked_list(graph, head, x0y1z_);
            extend_linked_list(graph, head, x0y_z0);
            extend_linked_list(graph, head, x0y_z1);
            extend_linked_list(graph, head, x1y0z_);
            extend_linked_list(graph, head, x1y1z_);
            extend_linked_list(graph, head, x1y_z0);
            extend_linked_list(graph, head, x1y_z1);
            extend_linked_list(graph, head, x_y0z0);
            extend_linked_list(graph, head, x_y0z1);
            extend_linked_list(graph, head, x_y1z0);
            extend_linked_list(graph, head, x_y1z1);
            extend_linked_list(graph, head, x0y0z0);
            extend_linked_list(graph, head, x0y0z1);
            extend_linked_list(graph, head, x0y1z0);
            extend_linked_list(graph, head, x0y1z1);
            extend_linked_list(graph, head, x1y0z0);
            extend_linked_list(graph, head, x1y0z1);
            extend_linked_list(graph, head, x1y1z0);
            extend_linked_list(graph, head, x1y1z1);

            head = pop_linked_list(head);
        }
    }

    status = write_bin(path_out, &clusters);
    if (status != 0) { return status; }

    free(smif.data);
    free(clusters.data);
    free(graph);

    printf(
        ">>> Wrote %i clusters to %s.\n... resolution: (nx=%u ny=%u nz=%u)\n... threshold: %g\n",
        (int)current_cluster, path_out, clusters.nx, clusters.ny, clusters.nz, (double)thr
    );
    return 0;
}

// -----------------------------------------------------------------------------
