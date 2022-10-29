/*
 * The author of this software is Steven Fortune.  Copyright (c) 1994 by AT&T
 * Bell Laboratories.
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

/* 
 * This code was originally written by Stephan Fortune in C code.  I, Shane O'Sullivan,
 * have since modified it, encapsulating it in a C++ class and, fixing memory leaks and
 * adding accessors to the Voronoi Edges.
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

/* 
 * Java Version by Zhenyu Pan
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

/*
 * New functions have been added by Quatja
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

/*
 * Modifications have been made and new functions have been added by Anlan Zhang
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

/*
 * CPP Version by Anlan Zhang
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

// Created by Anlan Zhang, 10/26/2022
#ifndef __VORONOI_HPP__
#define __VORONOI_HPP__

#include<vector>

class Point {
public:
    double x, y;

    Point();
    Point(double x, double y);
    ~Point();
    void set_point(double x, double y);
    bool equal_to(Point p);
};

class Site {
public:
    Point coord;
    int site_nbr;

    Site();
    ~Site();
};

class Edge {
public:
    double a;
    double b;
    double c;
    Site *ep[2]; // End points
    Site *reg[2]; // Sites this edge bisects
    int edge_nbr;

    Edge();
    ~Edge();
};

class GraphEdge {
public:
    double x1, y1, x2, y2;
    int site1, site2;

    GraphEdge();
    ~GraphEdge();
};

class HalfEdge {
public:
    HalfEdge *el_left;
    HalfEdge *el_right;
    Edge *el_edge;
    bool deleted;
    int el_pm;
    Site *vertex;
    double y_star;
    HalfEdge *pq_next;

    HalfEdge();
    ~HalfEdge();
};

class Voronoi {
private:
    double border_min_x;
    double border_max_x;
    double border_min_y;
    double border_max_y;

    int site_idx;

    double x_min;
    double y_min;

    double delta_x;
    double delta_y;

    int n_vertices;
    int n_edges;
    int n_sites;

    std::vector<Site *> sites;
    Site *bottom_site;

    int sqrt_n_sites;

    double min_distance_between_sites;

    int pq_count;
    int pq_min;
    int pq_hash_size;

    std::vector<HalfEdge *> pq_hash;

    const int LE = 0;
    const int RE = 1;

    int el_hash_size;
    std::vector<HalfEdge *> el_hash;
    HalfEdge *el_left_end, *el_right_end;
    std::vector<GraphEdge> all_edges;

    void add_corner_edges(const double x, const double y, const double max_y);

    void sort(double *x_values_in, double *y_values_in, int count);

    void qsort(std::vector<Site *> sites);

    void sort_node(double *x_values, double *y_values, int num_points);

    /* return a single in-storage site */
    Site *next_one();

    Edge *bisect(Site *s1, Site *s2);

    void make_vertex(Site *v);

    bool pq_initialize();

    int pq_bucket(HalfEdge *he);

    /* push the HalfEdge into the ordered linked list of vertices */
    void pq_insert(HalfEdge *he, Site *v, double offset);

    /* remove the HalfEdge from the list of vertices */
    void pq_delete(HalfEdge *he);

    bool pq_empty();

    Point get_pq_min();

    HalfEdge *pq_extract_min();

    HalfEdge *he_create(Edge *e, int pm);

    bool el_initialize();

    HalfEdge *el_right(HalfEdge *he);

    HalfEdge *el_left(HalfEdge *he);

    Site *left_reg(HalfEdge *he);

    void el_insert(HalfEdge *lb, HalfEdge *new_he);

    /*
     * This delete routine can't reclaim node, since pointers from hash table
     * may be present.
     */
    void el_delete(HalfEdge *he);

    /* Get entry from hash table, pruning any deleted nodes */
    HalfEdge *el_get_hash(int b);

    HalfEdge *el_left_bnd(Point p);

    void push_graph_edge(Site *left_site, Site *right_site, double x1, double y1, double x2, double y2);

    void clip_line(Edge *e);

    void endpoint(Edge *e, int lr, Site *s);

    /* returns 1 if p is to right of halfedge e */
    bool right_of(HalfEdge *el, Point p);

    Site *right_reg(HalfEdge *he);

    double dist(Site *s, Site *t);

    /* 
     * Create a new site where the HalfEdges el1 and el2 intersect - note that
     * the Point in the argument list is not used, don't know why it's there
     */
    Site *intersect(HalfEdge *el1, HalfEdge *el2);

    /*
     * implicit parameters: nsites, sqrt_nsites, xmin, xmax, ymin, ymax, deltax,
     * deltay (can all be estimates). Performance suffers if they are wrong;
     * better to make nsites, deltax, and deltay too big than too small. (?)
     */
    bool generate();
    
public:
    Voronoi(double min_distance_between_sites);
    ~Voronoi();

    /**
     * 
     * @param x_values_in Array of X values for each site.
     * @param y_values_in Array of Y values for each site. Must be identical length to yValuesIn
     * @param count size of x_values_in (or y_values_in)
     * @param min_x The minimum X of the bounding box around the voronoi
     * @param max_x The maximum X of the bounding box around the voronoi
     * @param min_y The minimum Y of the bounding box around the voronoi
     * @param max_y The maximum Y of the bounding box around the voronoi
     */
    void generate_voronoi(double *x_values_in, double *y_values_in, int count, double min_x, double max_x, double min_y, double max_y);

    // return the indexes of site of each pair of neighbors
    // noted that the index of site is equal to its coord's index in x_values_in/y_values_in
    std::vector<std::vector<int>> get_all_neighbor_pairs();

    void print_all_sites();
    void print_all_graph_edges();

};

#endif