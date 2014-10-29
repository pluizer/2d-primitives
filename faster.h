#ifndef __faster_h_
#define __faster_h_

extern void polygon_add(float* polygon, float* vect, int q_vertices, float* ret);

extern void polygon_sub(float* polygon, float* vect, int q_vertices, float* ret);

extern void polygon_rotate(float* polygon, float angle, float* origin, int q_vertices, float* ret);

#endif //__faster_h_
