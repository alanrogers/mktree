/** defined in alloc2d.c **/
void       **alloc2d(int dim1, int dim2, int size);
int          free2d(void **mat);
void      ***alloc3d(int dim1, int dim2, int dim3, int size);
int          free3d(void ***mat);
void     ****alloc4d(int dim1, int dim2, int dim3, int dim4, int size);
int          free4d(void ****mat);
void       **alloclt(int dim, int size);
void       **allocut(int dim, int size);
