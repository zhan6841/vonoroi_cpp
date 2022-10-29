#include "Vonoroi.hpp"

#include <stdio.h>
#include <algorithm>
#include <math.h>

// Point
Point::Point() {}

Point::~Point() {}

Point::Point(double x, double y) {
    this->x = x;
    this->y = y;
}

void Point::set_point(double x, double y) {
    this->x = x;
    this->y = y;
}

bool Point::equal_to(Point p) {
    return (this->x == p.x && this->y == p.y);
}

// Site
Site::Site() {}

Site::~Site() {}

// Edge
Edge::Edge() {
    this->ep[0] = NULL;
    this->ep[1] = NULL;
    this->reg[0] = NULL;
    this->reg[1] = NULL;
}

Edge::~Edge() {}

// GraphEdge
GraphEdge::GraphEdge() {}

GraphEdge::~GraphEdge() {}

// HalfEdge
HalfEdge::HalfEdge() {
    this->el_left = NULL;
    this->el_right = NULL;
    this->el_edge = NULL;
    this->vertex = NULL;
    this->pq_next = NULL;
    this->deleted = false;
}

HalfEdge::~HalfEdge() {};

// Comparator
static bool compare_site_less(Site *p1, Site *p2) {
    Point s1 = p1->coord, s2 = p2->coord;

    if (s1.y < s2.y) {
        return true;
    }

    if (s1.y > s2.y) {
        return false;
    }

    if (s1.x < s2.x) {
        return true;
    }

    if (s1.x > s2.x) {
        return false;
    }

    return true;
}

// Voronoi
/*********************************************************
 * Public methods
 ********************************************************/

Voronoi::Voronoi(double min_distance_between_sites) {
    // siteidx = 0;
    // sites = null;

    // allEdges = null;
    // this.minDistanceBetweenSites = minDistanceBetweenSites;

    this->site_idx = 0;
    
    this->min_distance_between_sites = min_distance_between_sites;

    this->bottom_site = NULL;
    this->el_left_end = NULL;
    this->el_right_end;
}

Voronoi::~Voronoi() {
    for (size_t i = 0; i < this->sites.size(); i++) {
        if (this->sites[i] != NULL) {
            delete this->sites[i];
        }
    }

    for (size_t i = 0; i < this->pq_hash.size(); i++) {
        if (this->pq_hash[i] != NULL) {
            delete this->pq_hash[i];
        }
    }

    for (size_t i = 0; i < this->el_hash.size(); i++) {
        if (this->el_hash[i] != NULL) {
            if (this->el_hash[i]->deleted) {
                delete this->el_hash[i];
            }
        }
    }

    if (this->el_left_end != NULL) {
        delete this->el_left_end;
    }

    if (this->el_right_end != NULL) {
        delete this->el_right_end;
    }
}

/**
 * 
 * @param x_values_in Array of X values for each site.
 * @param y_values_n Array of Y values for each site. Must be identical length to yValuesIn
 * @param count size of x_values_in (or y_values_in)
 * @param min_x The minimum X of the bounding box around the voronoi
 * @param max_x The maximum X of the bounding box around the voronoi
 * @param min_y The minimum Y of the bounding box around the voronoi
 * @param max_y The maximum Y of the bounding box around the voronoi
 */
void Voronoi::generate_voronoi(
        double *x_values_in, 
        double *y_values_in, 
        int count,
        double min_x, 
        double max_x, 
        double min_y, 
        double max_y) {
    if (x_values_in == NULL || y_values_in == NULL || count == 0) {
        return;
    }

    this->sort(x_values_in, y_values_in, count);

    this->border_min_x = min_x;
    this->border_min_y = min_y;
    this->border_max_x = max_x;
    this->border_max_y = max_y;

    this->site_idx = 0;
    this->generate();
}

// return the indexes of site of each pair of neighbors
// noted that the index of site is equal to its coord's index in x_values_in/y_values_in
std::vector<std::vector<int>> Voronoi::get_all_neighbor_pairs() {
    std::vector<std::vector<int>> all_neighbor_pairs;
    for (GraphEdge edge : this->all_edges) {
        std::vector<int> pair = {edge.site1, edge.site2};
        all_neighbor_pairs.push_back(pair);
    }
    return all_neighbor_pairs;
}

void Voronoi::print_all_sites() {
    for (Site *site : this->sites) {
        fprintf(stdout, "%d: %f, %f\n", site->site_nbr, site->coord.x, site->coord.y);
    }
}

void Voronoi::print_all_graph_edges() {
    for (GraphEdge edge : this->all_edges) {
        fprintf(stdout, "%f, %f; %f, %f; %d, %d\n", edge.x1, edge.y1, edge.x2, edge.y2, edge.site1, edge.site2);
    }
}

/*********************************************************
 * Private methods - implementation details
 ********************************************************/

void Voronoi::add_corner_edges(const double x, const double y, const double max_y) {
    double cur_y = y;
    double min_dist = INFINITY;

    // make a vertical line from the start corner to the closest region or the end corner
    for (size_t i = 0; i < this->all_edges.size(); i++) {
        GraphEdge edge = this->all_edges[i];

        if (edge.x1 == x) {
            double dist = std::abs(edge.y1 - y);
            if (dist < min_dist) {
                min_dist = dist;
                cur_y = edge.y1;
            }
        }
        if (edge.x2 == x) {
            double dist = std::abs(edge.y2 - y);
            if (dist < min_dist) {
                min_dist = dist;
                cur_y = edge.y2;
            }
        }
    }
    if (cur_y == y) {
        cur_y = max_y;
    }


    // find the closest site from the center of the line
    GraphEdge edge;
    edge.x1 = x;
    edge.y1 = y;
    edge.x2 = x;
    edge.y2 = cur_y;

    min_dist = INFINITY;
    int index = 0;
    double y_center = (edge.y2+edge.y1)/2;
    for (Site *site : sites) {
        double my_dist = std::sqrt(std::pow(site->coord.x - x, 2) + std::pow(site->coord.y - y_center, 2));
        if (my_dist < min_dist) {
            min_dist = my_dist;
            index = site->site_nbr;
        }
    }

    edge.site1 = index;
    edge.site2 = index;

    this->all_edges.push_back(edge);
}

void Voronoi::sort(double *x_values_in, double *y_values_in, int count) {
    // sites = null;
    // allEdges = new LinkedList<>();

    this->n_sites = count;
    this->n_vertices = 0;
    this->n_edges = 0;

    double sn = (double) this->n_sites + 4;
    this->sqrt_n_sites = (int) sqrt(sn);

    // // Copy the inputs so we don't modify the originals
    // double *x_values = new double[count];
    // double *y_values = new double[count];
    // for (int i = 0; i < count; i++) {
    //     x_values[i] = x_values_in[i];
    //     y_values[i] = y_values_in[i];
    // }
    // this->sort_node(x_values, y_values, count);
    this->sort_node(x_values_in, y_values_in, count);
}

void Voronoi::qsort(std::vector<Site *> sites) {
    std::sort(this->sites.begin(), this->sites.end(), compare_site_less);
}

void Voronoi::sort_node(double *x_values, double *y_values, int num_points) {
    int i;

    this->n_sites = num_points;
    this->sites = std::vector<Site *> (this->n_sites, NULL);
    this->x_min = x_values[0];
    this->y_min = y_values[0];
    
    double x_max = x_values[0];
    double y_max = y_values[0];
    for (i = 0; i < this->n_sites; i++) {
        sites[i] = new Site();
        sites[i]->coord.set_point(x_values[i], y_values[i]);
        sites[i]->site_nbr = i;

        if (x_values[i] < this->x_min) {
            this->x_min = x_values[i];
        } else if (x_values[i] > x_max) {
            x_max = x_values[i];
        }

        if (y_values[i] < this->y_min) {
            this->y_min = y_values[i];
        } else if (y_values[i] > y_max) {
            y_max = y_values[i];
        }
    }
    this->qsort(this->sites);
    this->delta_y = y_max - this->y_min;
    this->delta_x = x_max - this->x_min;
}

/* return a single in-storage site */
Site* Voronoi::next_one() {
    Site *s = NULL;
    
    if (this->site_idx < this->n_sites) {
        s = this->sites[this->site_idx];
        this->site_idx += 1;
        return (s);
    } else {
        return (NULL);
    }
}

Edge* Voronoi::bisect(Site *s1, Site *s2) {
    double dx;
    double dy;
    double adx;
    double ady;
    Edge *newedge = NULL;

    newedge = new Edge();

    // store the sites that this edge is bisecting
    newedge->reg[0] = s1;
    newedge->reg[1] = s2;
    // to begin with, there are no endpoints on the bisector - it goes to infinity
    newedge->ep[0] = NULL;
    newedge->ep[1] = NULL;

    // get the difference in x dist between the sites
    dx = s2->coord.x - s1->coord.x;
    dy = s2->coord.y - s1->coord.y;
    // make sure that the difference in positive
    adx = dx > 0 ? dx : -dx;
    ady = dy > 0 ? dy : -dy;
    newedge->c = (double) (s1->coord.x * dx + s1->coord.y * dy + (dx * dx + dy * dy) * 0.5);// get the slope of the line

    if (adx > ady) {
        newedge->a = 1.0f;
        newedge->b = dy / dx;
        newedge->c /= dx;// set formula of line, with x fixed to 1
    } else {
        newedge->b = 1.0f;
        newedge->a = dx / dy;
        newedge->c /= dy;// set formula of line, with y fixed to 1
    }

    newedge->edge_nbr = this->n_edges;

    this->n_edges += 1;
    return (newedge);
}

void Voronoi::make_vertex(Site *v) {
    v->site_nbr = this->n_vertices;
    this->n_vertices += 1;
}

bool Voronoi::pq_initialize() {
    this->pq_count = 0;
    this->pq_min = 0;
    this->pq_hash_size = 4 * this->sqrt_n_sites;
    this->pq_hash = std::vector<HalfEdge *> (this->pq_hash_size, NULL);

    for (int i = 0; i < this->pq_hash_size; i += 1) {
        this->pq_hash[i] = new HalfEdge();
    }
    return true;
}

int Voronoi::pq_bucket(HalfEdge *he) {
    int bucket;

    bucket = (int) ((he->y_star - this->y_min) / this->delta_y * this->pq_hash_size);
    if (bucket < 0) {
        bucket = 0;
    }
    if (bucket >= this->pq_hash_size) {
        bucket = this->pq_hash_size - 1;
    }
    if (bucket < this->pq_min) {
        this->pq_min = bucket;
    }
    return (bucket);
}

/* push the HalfEdge into the ordered linked list of vertices */
void Voronoi::pq_insert(HalfEdge *he, Site *v, double offset) {
    HalfEdge *last = NULL, *next = NULL;

    he->vertex = v;
    he->y_star = (double) (v->coord.y + offset);
    last = this->pq_hash[this->pq_bucket(he)];
    while ((next = last->pq_next) != NULL
            && (he->y_star > next->y_star || 
            (he->y_star == next->y_star && v->coord.x > next->vertex->coord.x))) {
        last = next;
    }
    he->pq_next = last->pq_next;
    last->pq_next = he;
    this->pq_count += 1;
}

/* remove the HalfEdge from the list of vertices */
void Voronoi::pq_delete(HalfEdge *he) {
    HalfEdge *last = NULL;

    if (he->vertex != NULL) {
        last = this->pq_hash[this->pq_bucket(he)];
        while (last->pq_next != he) {
            last = last->pq_next;
        }

        last->pq_next = he->pq_next;
        this->pq_count -= 1;
        he->vertex = NULL;
        delete he;
    }
}

bool Voronoi::pq_empty() {
    return (this->pq_count == 0);
}

Point Voronoi::get_pq_min() {
    Point answer;

    while (this->pq_hash[this->pq_min]->pq_next == NULL) {
        this->pq_min += 1;
    }
    answer.x = this->pq_hash[this->pq_min]->pq_next->vertex->coord.x;
    answer.y = this->pq_hash[this->pq_min]->pq_next->y_star;
    return (answer);
}

HalfEdge* Voronoi::pq_extract_min() {
    HalfEdge *curr = NULL;

    curr = this->pq_hash[this->pq_min]->pq_next;
    this->pq_hash[this->pq_min]->pq_next = curr->pq_next;
    this->pq_count -= 1;
    return (curr);
}

HalfEdge* Voronoi::he_create(Edge *e, int pm) {
    HalfEdge *answer = NULL;
    answer = new HalfEdge();
    answer->el_edge = e;
    answer->el_pm = pm;
    answer->pq_next = NULL;
    answer->vertex = NULL;
    return (answer);
}

bool Voronoi::el_initialize() {
    int i;
    this->el_hash_size = 2 * this->sqrt_n_sites;
    this->el_hash = std::vector<HalfEdge *> (this->el_hash_size, NULL);

    for (i = 0; i < this->el_hash_size; i += 1) {
        this->el_hash[i] = NULL;
    }
    this->el_left_end = this->he_create(NULL, 0);
    this->el_right_end = this->he_create(NULL, 0);
    this->el_left_end->el_left = NULL;
    this->el_left_end->el_right = this->el_right_end;
    this->el_right_end->el_left = this->el_left_end;
    this->el_right_end->el_right = NULL;
    this->el_hash[0] = this->el_left_end;
    this->el_hash[this->el_hash_size - 1] = this->el_right_end;

    return true;
}

HalfEdge* Voronoi::el_right(HalfEdge *he) {
    return (he->el_right);
}

HalfEdge* Voronoi::el_left(HalfEdge *he) {
    return (he->el_left);
}

Site* Voronoi::left_reg(HalfEdge *he) {
    if (he->el_edge == NULL) {
        return (this->bottom_site);
    }
    return (he->el_pm == this->LE ? he->el_edge->reg[this->LE] : he->el_edge->reg[this->RE]);
}

void Voronoi::el_insert(HalfEdge *lb, HalfEdge *new_he) {
    new_he->el_left = lb;
    new_he->el_right = lb->el_right;
    (lb->el_right)->el_left = new_he;
    lb->el_right = new_he;
}

/*
 * This delete routine can't reclaim node, since pointers from hash table
 * may be present.
 */
void Voronoi::el_delete(HalfEdge *he) {
    (he->el_left)->el_right = he->el_right;
    (he->el_right)->el_left = he->el_left;
    he->deleted = true;
}

/* Get entry from hash table, pruning any deleted nodes */
HalfEdge* Voronoi::el_get_hash(int b) {
    HalfEdge *he = NULL;

    if (b < 0 || b >= this->el_hash_size) {
        return (NULL);
    }
    he = this->el_hash[b];
    if (he == NULL || !he->deleted) {
        return (he);
    }

    /* Hash table points to deleted half edge. Patch as necessary. */
    this->el_hash[b] = NULL;
    return (NULL);
}

HalfEdge* Voronoi::el_left_bnd(Point p) {
    int i, bucket;
    HalfEdge *he = NULL;

    /* Use hash table to get close to desired halfedge */
    // use the hash function to find the place in the hash map that this
    // HalfEdge should be
    bucket = (int) ((p.x - this->x_min) / this->delta_x * this->el_hash_size);

    // make sure that the bucket position in within the range of the hash array
    if (bucket < 0) {
        bucket = 0;
    }
    if (bucket >= this->el_hash_size) {
        bucket = this->el_hash_size - 1;
    }

    he = this->el_get_hash(bucket);
    // if the HE isn't found, search backwards and forwards in the hash map
    // for the first non-null entry
    if (he == NULL) {
        for (i = 1; i < this->el_hash_size; i += 1) {
            if ((he = this->el_get_hash(bucket - i)) != NULL) {
                break;
            }
            if ((he = this->el_get_hash(bucket + i)) != NULL) {
                break;
            }
        }
    }
    /* Now search linear list of halfedges for the correct one */
    if (he == this->el_left_end || (he != this->el_right_end && this->right_of(he, p))) {
        // keep going right on the list until either the end is reached, or
        // you find the 1st edge which the point isn't to the right of
        do {
            he = he->el_right;
        } while (he != this->el_right_end && this->right_of(he, p));
        he = he->el_left;
    } else {
        // if the point is to the left of the HalfEdge, then search left for
        // the HE just to the left of the point
        do {
            he = he->el_left;
        } while (he != this->el_left_end && !this->right_of(he, p));
    }

    /* Update hash table and reference counts */
    if (bucket > 0 && bucket < this->el_hash_size - 1) {
        this->el_hash[bucket] = he;
    }
    return (he);
}

void Voronoi::push_graph_edge(Site *left_site, Site *right_site, double x1, double y1, double x2, double y2) {
    GraphEdge new_edge;
    new_edge.x1 = x1;
    new_edge.y1 = y1;
    new_edge.x2 = x2;
    new_edge.y2 = y2;

    new_edge.site1 = left_site->site_nbr;
    new_edge.site2 = right_site->site_nbr;
    this->all_edges.push_back(new_edge);
}

void Voronoi::clip_line(Edge *e) {
    Site *s1 = NULL;
    Site *s2 = NULL;

    double pxmin;
    double pxmax;
    double pymin;
    double pymax;

    double x1;
    double x2;
    double y1;
    double y2;

    x1 = e->reg[0]->coord.x;
    x2 = e->reg[1]->coord.x;
    y1 = e->reg[0]->coord.y;
    y2 = e->reg[1]->coord.y;

    // if the distance between the two points this line was created from is
    // less than the square root of 2, then ignore it
    if (sqrt(((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1))) < this->min_distance_between_sites) {
        return;
    }
    pxmin = this->border_min_x;
    pxmax = this->border_max_x;
    pymin = this->border_min_y;
    pymax = this->border_max_y;

    if (e->a == 1.0 && e->b >= 0.0) {
        s1 = e->ep[1];
        s2 = e->ep[0];
    } else {
        s1 = e->ep[0];
        s2 = e->ep[1];
    }

    if (e->a == 1.0) {
        y1 = pymin;
        if (s1 != NULL && s1->coord.y > pymin) {
            y1 = s1->coord.y;
        }
        if (y1 > pymax) {
            y1 = pymax;
        }
        x1 = e->c - e->b * y1;
        y2 = pymax;
        if (s2 != NULL && s2->coord.y < pymax) {
            y2 = s2->coord.y;
        }

        if (y2 < pymin) {
            y2 = pymin;
        }
        x2 = (e->c) - (e->b) * y2;
        if (((x1 > pxmax) & (x2 > pxmax)) | ((x1 < pxmin) & (x2 < pxmin))) {
            return;
        }
        if (x1 > pxmax) {
            x1 = pxmax;
            y1 = (e->c - x1) / e->b;
        }
        if (x1 < pxmin) {
            x1 = pxmin;
            y1 = (e->c - x1) / e->b;
        }
        if (x2 > pxmax) {
            x2 = pxmax;
            y2 = (e->c - x2) / e->b;
        }
        if (x2 < pxmin) {
            x2 = pxmin;
            y2 = (e->c - x2) / e->b;
        }
    } else {
        x1 = pxmin;
        if (s1 != NULL && s1->coord.x > pxmin) {
            x1 = s1->coord.x;
        }
        if (x1 > pxmax) {
            x1 = pxmax;
        }
        y1 = e->c - e->a * x1;
        x2 = pxmax;
        if (s2 != NULL && s2->coord.x < pxmax) {
            x2 = s2->coord.x;
        }
        if (x2 < pxmin) {
            x2 = pxmin;
        }
        y2 = e->c - e->a * x2;
        if (((y1 > pymax) & (y2 > pymax)) | ((y1 < pymin) & (y2 < pymin))) {
            return;
        }
        if (y1 > pymax) {
            y1 = pymax;
            x1 = (e->c - y1) / e->a;
        }
        if (y1 < pymin) {
            y1 = pymin;
            x1 = (e->c - y1) / e->a;
        }
        if (y2 > pymax) {
            y2 = pymax;
            x2 = (e->c - y2) / e->a;
        }
        if (y2 < pymin) {
            y2 = pymin;
            x2 = (e->c - y2) / e->a;
        }
    }

    this->push_graph_edge(e->reg[0], e->reg[1], x1, y1, x2, y2);
}

void Voronoi::endpoint(Edge *e, int lr, Site *s) {
    e->ep[lr] = s;
    if (e->ep[this->RE - lr] == NULL){
        return;
    }
    this->clip_line(e);
}

/* returns 1 if p is to right of halfedge e */
bool Voronoi::right_of(HalfEdge *el, Point p) {
    Edge *e = NULL;
    Site *top_site = NULL;

    bool right_of_site;
    bool above;
    bool fast;

    double dxp;
    double dyp;
    double dxs;
    double t1;
    double t2;
    double t3;
    double yl;

    e = el->el_edge;
    top_site = e->reg[1];
    if (p.x > top_site->coord.x) {
        right_of_site = true;
    } else {
        right_of_site = false;
    }
        
    if (right_of_site && el->el_pm == this->LE) {
        return (true);
    }
    if (!right_of_site && el->el_pm == this->RE) {
        return (false);
    }

    if (e->a == 1.0) {
        dyp = p.y - top_site->coord.y;
        dxp = p.x - top_site->coord.x;
        fast = false;
        if ((!right_of_site & (e->b < 0.0)) | (right_of_site & (e->b >= 0.0))) {
            above = dyp >= e->b * dxp;
            fast = above;
        } else {
            above = (p.x + p.y * e->b > e->c);
            if (e->b < 0.0) {
                above = !above;
            }
            if (!above) {
                fast = true;
            }
        }
        if (!fast) {
            dxs = top_site->coord.x - (e->reg[0])->coord.x;
            above = (e->b * (dxp * dxp - dyp * dyp) < dxs * dyp
                    * (1.0 + 2.0 * dxp / dxs + e->b * e->b));
            if (e->b < 0.0) {
                above = !above;
            }
        }
    } else {
        /* e.b==1.0 */
        yl = e->c - e->a * p.x;
        t1 = p.y - yl;
        t2 = p.x - top_site->coord.x;
        t3 = yl - top_site->coord.y;
        above = (t1 * t1 > t2 * t2 + t3 * t3);
    }
    return (el->el_pm == this->LE ? above : !above);
}

Site* Voronoi::right_reg(HalfEdge *he) {
    // if this halfedge has no edge, return the bottom site (whatever that is)
    if (he->el_edge == NULL) {
        return (this->bottom_site);
    }

    // if the ELpm field is zero, return the site 0 that this edge bisects,
    // otherwise return site number 1
    return (he->el_pm == this->LE ? he->el_edge->reg[this->RE] : he->el_edge->reg[this->LE]);
}

double Voronoi::dist(Site *s, Site *t) {
    double dx;
    double dy;
    dx = s->coord.x - t->coord.x;
    dy = s->coord.y - t->coord.y;
    return sqrt(dx * dx + dy * dy);
}

/* 
 * Create a new site where the HalfEdges el1 and el2 intersect - note that
 * the Point in the argument list is not used, don't know why it's there
 */
Site* Voronoi::intersect(HalfEdge *el1, HalfEdge *el2) {
    Edge *e1 = NULL;
    Edge *e2 = NULL;
    Edge *e = NULL;

    HalfEdge *el = NULL;

    double d;
    double xint;
    double yint;

    bool right_of_site;

    Site *v = NULL;

    e1 = el1->el_edge;
    e2 = el2->el_edge;
    if (e1 == NULL || e2 == NULL) {
        return NULL;
    }

    // if the two edges bisect the same parent, return null
    if (e1->reg[1] == e2->reg[1]) {
        return NULL;
    }

    d = e1->a * e2->b - e1->b * e2->a;
    if (-1.0e-10 < d && d < 1.0e-10) {
        return NULL;
    }

    xint = (e1->c * e2->b - e2->c * e1->b) / d;
    yint = (e2->c * e1->a - e1->c * e2->a) / d;

    if ((e1->reg[1]->coord.y < e2->reg[1]->coord.y) || 
            (e1->reg[1]->coord.y == e2->reg[1]->coord.y && 
            e1->reg[1]->coord.x < e2->reg[1]->coord.x)) {
        el = el1;
        e = e1;
    } else {
        el = el2;
        e = e2;
    }

    right_of_site = xint >= e->reg[1]->coord.x;
    if ((right_of_site && el->el_pm == this->LE) || (!right_of_site && el->el_pm == this->RE)) {
        return NULL;
    }

    // create a new site at the point of intersection - this is a new vector
    // event waiting to happen
    v = new Site();
    v->coord.x = xint;
    v->coord.y = yint;
    return (v);
}

/*
 * implicit parameters: nsites, sqrt_nsites, xmin, xmax, ymin, ymax, deltax,
 * deltay (can all be estimates). Performance suffers if they are wrong;
 * better to make nsites, deltax, and deltay too big than too small. (?)
 */
bool Voronoi::generate() {
    Site *new_site = NULL;
    Site *bot = NULL;
    Site *top = NULL;
    Site *temp = NULL;
    Site *p = NULL;
    Site *v = NULL;

    Point new_int_star;

    int pm;

    HalfEdge *lbnd = NULL;
    HalfEdge *rbnd = NULL;
    HalfEdge *llbnd = NULL;
    HalfEdge *rrbnd = NULL;
    HalfEdge *bisector = NULL;

    Edge *e = NULL;

    this->pq_initialize();
    this->el_initialize();

    this->bottom_site = this->next_one();
    new_site = this->next_one();
    while (true) {
        if (!this->pq_empty()) {
            new_int_star = this->get_pq_min();
        }
        // if the lowest site has a smaller y value than the lowest vector
        // intersection,
        // process the site otherwise process the vector intersection

        if (new_site != NULL && 
                (this->pq_empty() || new_site->coord.y < new_int_star.y || 
                (new_site->coord.y == new_int_star.y && new_site->coord.x < new_int_star.x))) {
            /* new site is smallest -this is a site event */
            // get the first HalfEdge to the LEFT of the new site
            lbnd = this->el_left_bnd((new_site->coord));
            // get the first HalfEdge to the RIGHT of the new site
            rbnd = this->el_right(lbnd);
            // if this halfedge has no edge,bot =bottom site (whatever that is)
            bot = this->right_reg(lbnd);
            // create a new edge that bisects
            e = this->bisect(bot, new_site);


            // create a new HalfEdge, setting its ELpm field to 0
            bisector = this->he_create(e, LE);
            // insert this new bisector edge between the left and right
            // vectors in a linked list
            this->el_insert(lbnd, bisector);

            // if the new bisector intersects with the left edge,
            // remove the left edge's vertex, and put in the new one
            if ((p = this->intersect(lbnd, bisector)) != NULL) {
                this->pq_delete(lbnd);
                this->pq_insert(lbnd, p, dist(p, new_site));
            }
            lbnd = bisector;
            // create a new HalfEdge, setting its ELpm field to 1
            bisector = this->he_create(e, RE);
            // insert the new HE to the right of the original bisector
            // earlier in the IF stmt
            this->el_insert(lbnd, bisector);

            // if this new bisector intersects with the new HalfEdge
            if ((p = this->intersect(bisector, rbnd)) != NULL) {
                // push the HE into the ordered linked list of vertices
                this->pq_insert(bisector, p, dist(p, new_site));
            }
            new_site = this->next_one();
        } else if (!this->pq_empty()) {
            /* intersection is smallest - this is a vector event */
            // pop the HalfEdge with the lowest vector off the ordered list
            // of vectors
            lbnd = this->pq_extract_min();
            // get the HalfEdge to the left of the above HE
            llbnd = this->el_left(lbnd);
            // get the HalfEdge to the right of the above HE
            rbnd = this->el_right(lbnd);
            // get the HalfEdge to the right of the HE to the right of the
            // lowest HE
            rrbnd = this->el_right(rbnd);
            // get the Site to the left of the left HE which it bisects
            bot = this->left_reg(lbnd);
            // get the Site to the right of the right HE which it bisects
            top = this->right_reg(rbnd);

            v = lbnd->vertex; // get the vertex that caused this event
            this->make_vertex(v); // set the vertex number - couldn't do this
            // earlier since we didn't know when it would be processed
            this->endpoint(lbnd->el_edge, lbnd->el_pm, v);
            // set the endpoint of
            // the left HalfEdge to be this vector
            this->endpoint(rbnd->el_edge, rbnd->el_pm, v);
            // set the endpoint of the right HalfEdge to
            // be this vector
            this->el_delete(lbnd); // mark the lowest HE for
            // deletion - can't delete yet because there might be pointers
            // to it in Hash Map
            this->pq_delete(rbnd);
            // remove all vertex events to do with the right HE
            this->el_delete(rbnd); // mark the right HE for
            // deletion - can't delete yet because there might be pointers
            // to it in Hash Map
            pm = LE; // set the pm variable to zero

            if (bot->coord.y > top->coord.y) { 
                // if the site to the left of the event is higher than the Site
                // to the right of it, then swap them and set the 'pm' variable to 1
                temp = bot;
                bot = top;
                top = temp;
                pm = RE;
            }
            e = this->bisect(bot, top); // create an Edge (or line)
            // that is between the two Sites. This creates the formula of
            // the line, and assigns a line number to it
            bisector = this->he_create(e, pm); // create a HE from the Edge 'e',
            // and make it point to that edge
            // with its ELedge field
            this->el_insert(llbnd, bisector); // insert the new bisector to the
            // right of the left HE
            this->endpoint(e, RE - pm, v); // set one endpoint to the new edge
            // to be the vector point 'v'.
            // If the site to the left of this bisector is higher than the
            // right Site, then this endpoint
            // is put in position 0; otherwise in pos 1

            // if left HE and the new bisector intersect, then delete
            // the left HE, and reinsert it
            if ((p = this->intersect(llbnd, bisector)) != NULL) {
                this->pq_delete(llbnd);
                this->pq_insert(llbnd, p, dist(p, bot));
            }

            // if right HE and the new bisector intersect, then
            // reinsert it
            if ((p = this->intersect(bisector, rrbnd)) != NULL) {
                this->pq_insert(bisector, p, dist(p, bot));
            }
        } else {
            break;
        }
    }

    for (lbnd = this->el_right(this->el_left_end); lbnd != this->el_right_end; lbnd = this->el_right(lbnd)) {
        e = lbnd->el_edge;
        this->clip_line(e);
    }

    return true;
}