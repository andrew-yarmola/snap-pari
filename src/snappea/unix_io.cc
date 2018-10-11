/* 
   As you may have noticed, SnapPea 2.0 uses a new file format.
   It reads both the old and new formats, but writes only the new one.
   The purpose of this e-mail message is to provide you with some
   sample code so your programs can handle the new format, too.

   The protocol by which your program passes files to SnapPea's
   computational kernel has changed.  On the one hand I wanted
   to maintain SnapPea's convention of using text files to store
   its data (so people using different platforms can easily send
   files back and forth).  On the other hand, I didn't want the
   kernel to be hard-coded to require a pointer to a unix-style FILE.
   The solution is to let the UI (i.e. your program) read the text
   file in whatever manner it wants, and then pass the raw data to
   the kernel in a TriangulationData structure.  The kernel function
   data_to_triangulation() converts the raw TriangulationData to a
   full-fledged Triangulation.  [triangulation_io.h defines the
   TriangulationData structure and "TriangulationFileFormat"
   explains the file format itself, but you don't have to read either
   to use the sample code below.]

   [As an added bonus, if you are writing a program which creates
   its own manifolds on the fly, you no longer have to mess with the
   details of setting up a Triangulation structure.  You just
   fill out the minimal amount of information in the TriangulationData,
   and let the kernel do the grunt work.  (At present you needn't
   specify the hyperbolic structure or the peripheral curves unless
   you want, and it will be easy for me to set it up so that you needn't
   specify the cusps either.  Please let me know if you're in a hurry
   for that capability.)]
   
   Anyhow . . . getting on with the topic at hand, here's some sample
   code to read a manifold from an already open unix-style FILE.
   Note that ReadTriangulationFromFile() and ReadNewFileFormat()
   will be part of your program, while data_to_triangulation() is
   provided by the kernel, and read_old_manifold() is contained in
   OldFileFormat.hack.c.  (OldFileFormat.hack.c is part of the Mac UI
   code -- if you don't have that I'd be happy to e-mail it to you.)
   */

#include "unix_io.h"
#include "kernel.h"
#include <cstdlib>

using std::string;
using std::exit;

static FuncResult ReadNewFileFormat(FILE *fp, Triangulation *manifold);
static FuncResult read_old_manifold(FILE* fp, Triangulation* manifold);
static FuncResult read_the_file(FILE *fp, Triangulation *manifold);
static void fix_orientation(Triangulation *manifold);
static void reverse_meridians_where_necessary(Triangulation *manifold);
static void WriteNewFileFormat(FILE *fp, TriangulationData *data);

static short num_census_manifolds[8] = {0, 0, 0, 0, 0, 415, 962, 3552};

#define REG_TET_VOL 1.014941606409653625021202554

FILE* locate_file(const string& path, const string& name, char *mode)
{
  string dirs[10], path_name; 

  int i=0,j,n=0;

  /* Split path into separate words. */ 
  while (n<10) {
    while (i < path.size() && path[i]==' ') ++i; /* Skip spaces. */
    if (i==path.size()) break; 
    j = i; 
    while (i < path.size() && path[i]!=' ') ++i; /* Go to the end of the word. */
    dirs[n] = path.substr(j,i-j); 
    ++n; 
  }

  FILE* fp; 
  for (j=0; j<n; j++) {
    path_name = dirs[j] + "/" + name;
    fp = fopen(path_name.c_str(), mode);
    if (fp) return fp; 
  }
  return NULL; 
}

#if 0
Triangulation* read_manifold_file(const char* path, const char* filename)
{
  Triangulation* manifold = NEW_STRUCT(Triangulation); 
  if (!manifold) return NULL;

  if (read_manifold_file(path, filename, manifold)==func_OK)
    return manifold;

  my_free(manifold);
  return NULL;
}

// Uses the name of the file as the name of the manifold if it
// is an old format file. 

FuncResult read_manifold_file(
    const char      *path, 
    const char      *filename,
    Triangulation   *manifold)
{
  Boolean         newFormat;
  FuncResult      result; 

  FILE *fp = locate_file(path, filename, "r"); 
  if (fp == NULL) {
    printf("File %s not found\n", (char*)filename);
    return func_failed;
  }

  //  Take a peek at the first line to see whether this is
  //  the new file format or the old one.
  newFormat = (getc(fp) == '%');
  rewind(fp);

  if (newFormat) {
    result = ReadNewFileFormat(fp, manifold);
    fclose(fp); 
    return result;
  }

  // else old file format
  result = read_old_manifold(fp, manifold);
  fclose(fp);

  if (result == func_OK)
    set_triangulation_name(manifold, (char*)filename);

  return result;
}
#endif

static FuncResult read_manifold_file(
    FILE            *fp,
    Triangulation   *manifold)
{
    Boolean         newFormat;

    //  Take a peek at the first line to see whether this is
    //  the new file format or the old one.
    newFormat = (getc(fp) == '%');
    rewind(fp);

    if (newFormat)
	return ReadNewFileFormat(fp, manifold);
    // else
    return read_old_manifold(fp, manifold);
}

Triangulation* read_manifold_file(FILE* fp)
{
  Triangulation* manifold = NEW_STRUCT(Triangulation); 
  if (!manifold) return NULL;

  if (read_manifold_file(fp, manifold)==func_OK)
    return manifold;

  my_free(manifold);
  return NULL;
}


/* 
   Here's an implementation of ReadNewFileFormat().
   In spite of its length, the basic idea is simple:  it just reads
   through the FILE copying the data to a TriangulationData structure.
*/

static FuncResult ReadNewFileFormat(FILE *fp, Triangulation *manifold)
{
    TriangulationData   theData;
    char                theScratchString[100];
    int                 theTotalNumCusps,
			i,
			j,
			k,
			v,
			f;

    //  Read and ignore the header (% Triangulation).
    fgets(theScratchString, 100, fp);

    //  Initialize the TriangulationData.
    theData.name              = NULL;
    theData.cusp_data         = NULL;
    theData.tetrahedron_data  = NULL;
    
    //  Allocate and read the manifold's name.
    theData.name = new char [100];
    if (!theData.name) goto error;

    //  The name will be on the first nonempty line.
    do
      fgets(theData.name, 100, fp);
    while (theData.name[0] == '\n');
    //  Overwrite the newline character.
    theData.name[strlen(theData.name) - 1] = 0;
    
    //  Read the filled solution type.
    fscanf(fp, "%s", theScratchString);
    if (strcmp(theScratchString, "not_attempted") == 0)
      theData.solution_type = not_attempted;
    else if (strcmp(theScratchString, "geometric_solution") == 0)
      theData.solution_type = geometric_solution;
    else if (strcmp(theScratchString, "nongeometric_solution") == 0)
      theData.solution_type = nongeometric_solution;
    else if (strcmp(theScratchString, "flat_solution") == 0)
      theData.solution_type = flat_solution;
    else if (strcmp(theScratchString, "degenerate_solution") == 0)
      theData.solution_type = degenerate_solution;
    else if (strcmp(theScratchString, "other_solution") == 0)
      theData.solution_type = other_solution;
    else if (strcmp(theScratchString, "no_solution") == 0)
      theData.solution_type = no_solution;
    else {
      fprintf(stderr, "solution type not recognized\n");
      goto error;
    }

    //  Read the volume.
    fscanf(fp, "%lf", &theData.volume);
    
    //  Read the orientability.
    fscanf(fp, "%s", theScratchString);
    if (strcmp(theScratchString, "oriented_manifold") == 0)
      theData.orientability = oriented_manifold;
    else if (strcmp(theScratchString, "nonorientable_manifold") == 0)
      theData.orientability = nonorientable_manifold;
    else {
      fprintf(stderr, "orientability not recognized\n");
      goto error;
    }
    
    //  Read the Chern-Simons invariant, if present.
    fscanf(fp, "%s", theScratchString);
    if (strcmp(theScratchString, "CS_known") == 0)
      theData.CS_value_is_known = TRUE;
    else if (strcmp(theScratchString, "CS_unknown") == 0)
      theData.CS_value_is_known = FALSE;
    else {
      fprintf(stderr, "Chern-Simons entry not recognized\n");
      goto error;
    }

    if (theData.CS_value_is_known == TRUE)
      fscanf(fp, "%lf", &theData.CS_value);
    else
      theData.CS_value = 0.0;
    
    //  Read the number of cusps, allocate an array for the cusp data,
    //  and read the cusp data.
    fscanf(fp, "%d%d",
	   &theData.num_or_cusps,
	   &theData.num_nonor_cusps);
    theTotalNumCusps = theData.num_or_cusps
      + theData.num_nonor_cusps;
    theData.cusp_data = new CuspData [theTotalNumCusps];
    if (!theData.cusp_data) goto error;
    for (i = 0; i < theTotalNumCusps; i++)
      {
	if (fscanf(fp, "%s%lf%lf",
		   theScratchString,
		   &theData.cusp_data[i].m,
		   &theData.cusp_data[i].l) != 3) {
	  fprintf(stderr, "cusp data not readable\n");
	  goto error;
	}
	switch (theScratchString[0])
	  {
	  case 't':
	  case 'T':
	    theData.cusp_data[i].topology = torus_cusp;
	    break;
	    
	  case 'k':
	  case 'K':
	    theData.cusp_data[i].topology = Klein_cusp;
	    break;
	    
	  default:
	    fprintf(stderr, "cusp topology not recognized\n");
	    goto error;
	  }
      }
    
    //  Read the number of tetrahedra, allocate an array for the
    //  tetrahedron data, and read the tetrahedron data.
    fscanf(fp, "%d", &theData.num_tetrahedra);
    theData.tetrahedron_data = new TetrahedronData [theData.num_tetrahedra];
    if (!theData.tetrahedron_data) goto error;
    for (i = 0; i < theData.num_tetrahedra; i++)
      {
	//  Read the neighbor indices.
	for (j = 0; j < 4; j++)
	  {
	    fscanf(fp, "%d",&theData.tetrahedron_data[i].neighbor_index[j]);
	    if (theData.tetrahedron_data[i].neighbor_index[j] < 0
		|| theData.tetrahedron_data[i].neighbor_index[j] >=
		theData.num_tetrahedra) {
	      fprintf(stderr, "invalid tetrahedron neighbor index found\n");
	      goto error;
	    }
	  }
	
	//  Read the gluings.
	for (j = 0; j < 4; j++)
	  for (k = 0; k < 4; k++)
	    {
	      fscanf(fp, "%1d",
		     &theData.tetrahedron_data[i].gluing[j][k]);
	      if
		(theData.tetrahedron_data[i].gluing[j][k] < 0 ||
		 theData.tetrahedron_data[i].gluing[j][k] > 3) {
		  fprintf(stderr, "invalid gluing data found\n");
		  goto error;
		}
	    }
	
	//  Read the cusp indices.
	for (j = 0; j < 4; j++)
	  {
	    fscanf(fp, "%d", &theData.tetrahedron_data[i].cusp_index[j]);
	    if (theData.tetrahedron_data[i].cusp_index[j] < 0
		|| theData.tetrahedron_data[i].cusp_index[j] >= theTotalNumCusps) {
	      fprintf(stderr, "invalid cusp index found\n");
	      goto error;
	    }
	  }
	
	//  Read the peripheral curves.
	for (j = 0; j < 2; j++)         //  meridian, longitude
	  for (k = 0; k < 2; k++)     //  righthanded, lefthanded
	    for (v = 0; v < 4; v++)
	      for (f = 0; f < 4; f++)
		fscanf(fp, "%d", &theData.tetrahedron_data[i].curve[j][k][v][f]);
	
	//  Read the filled shape (which the kernel ignores).
	fscanf(fp, "%lf%lf",
	       &theData.tetrahedron_data[i].filled_shape.real,
	       &theData.tetrahedron_data[i].filled_shape.imag);
      }

    triangulation_from_data(&theData, manifold);

    return func_OK;

  error:

    // these are OK on 0 pointers. 
    delete [] theData.name;
    delete [] theData.cusp_data;
    delete [] theData.tetrahedron_data;

    return func_failed;
}

/* **********************************************************************
   
   To write a Triangulation to a file, we just reverse the above process.
   (Note that when we allocate the TriangulationData, we must free it.
   But when the kernel allocates it, we must pass it back to the kernel
   function free_triangulation_data() to free it.  Otherwise the kernel
   will complain about memory leaks.)
   */

void write_manifold_file(
    FILE            *fp,
    Triangulation   *manifold)
{
    TriangulationData   *theTriangulationData;
  
    triangulation_to_data(manifold, &theTriangulationData);
    WriteNewFileFormat(fp, theTriangulationData);
    free_triangulation_data(theTriangulationData);
}


static void WriteNewFileFormat(
    FILE                *fp,
    TriangulationData   *data)
{
    int i,
	j,
	k,
	v,
	f;

    fprintf(fp, "%% Triangulation\n");

    if (data->name != NULL)
	fprintf(fp, "%s\n", data->name);
    else
	fprintf(fp, "untitled");

    switch (data->solution_type)
    {
	case not_attempted:
	    fprintf(fp, "not_attempted");
	    break;

	case geometric_solution:
	    fprintf(fp, "geometric_solution");
	    break;

	case nongeometric_solution:
	    fprintf(fp, "nongeometric_solution");
	    break;

	case flat_solution:
	    fprintf(fp, "flat_solution");
	    break;

	case degenerate_solution:
	    fprintf(fp, "degenerate_solution");
	    break;

	case other_solution:
	    fprintf(fp, "other_solution");
	    break;

	case no_solution:
	    fprintf(fp, "no_solution");
	    break;
    }

    fprintf(fp, "  %.8lf\n", data->volume);

    switch (data->orientability)
    {
	case oriented_manifold:
	    fprintf(fp, "oriented_manifold\n");
	    break;

	case nonorientable_manifold:
	    fprintf(fp, "nonorientable_manifold\n");
	    break;
    }

    if (data->CS_value_is_known == TRUE)
	fprintf(fp, "CS_known %.16lf\n", data->CS_value);
    else
	fprintf(fp, "CS_unknown\n");

    fprintf(fp, "\n%d %d\n", data->num_or_cusps, data->num_nonor_cusps);
    for (i = 0; i < data->num_or_cusps + data->num_nonor_cusps; i++)
	fprintf(fp, "    %s %16.12lf %16.12lf\n",
	    (data->cusp_data[i].topology == torus_cusp) ? "torus" : "Klein",
	    data->cusp_data[i].m,
	    data->cusp_data[i].l);
    fprintf(fp, "\n");

    fprintf(fp, "%d\n", data->num_tetrahedra);
    for (i = 0; i < data->num_tetrahedra; i++)
    {
	for (j = 0; j < 4; j++)
	    fprintf(fp, "%4d ", data->tetrahedron_data[i].neighbor_index[j]);
	fprintf(fp, "\n");

	for (j = 0; j < 4; j++)
	{
	    fprintf(fp, " ");
	    for (k = 0; k < 4; k++)
		fprintf(fp, "%d", data->tetrahedron_data[i].gluing[j][k]);
	}
	fprintf(fp, "\n");

	for (j = 0; j < 4; j++)
	    fprintf(fp, "%4d ", data->tetrahedron_data[i].cusp_index[j]);
	fprintf(fp, "\n");

	for (j = 0; j < 2; j++)         //  meridian, longitude
	    for (k = 0; k < 2; k++)     //  righthanded, lefthanded
	    {
		for (v = 0; v < 4; v++)
		    for (f = 0; f < 4; f++)
			fprintf(fp, " %2d", data->tetrahedron_data[i].curve[j][k][v][f]);
		fprintf(fp, "\n");
	    }

	fprintf(fp, "%16.12lf %16.12lf\n\n",
	    data->tetrahedron_data[i].filled_shape.real,
	    data->tetrahedron_data[i].filled_shape.imag);
    }
}

/*
   This is a hacked together file to read
   a manifold in the old file format and put it into the new Triangulation
   data structure.  The sole purpose of this file is to allow testing
   of the new program before it's capable of creating its own manifolds
   and/or converting old ones.  It cuts a lot of corners.  For example,
   it assumes the solution has maximal precision, it assumes the arguments
   of the TetShapes are in the expected range, it assumes the structure
   is a geometric_solution, etc.  To correct these deficiencies, it
   calls find_complete_hyperbolic_structure() to polish up whatever
   structure it started with.  This also hides the fact that the
   old snap pea didn't always provide very accurate solutions.

   This code should now handle orientations correctly, assuming the
   peripheral curves adhered to the standard orientation conventions
   to begin with.  E.g. link complements should all be OK.
*/

static FuncResult read_old_manifold(FILE* fp, Triangulation* manifold)
{

  initialize_triangulation(manifold);   /* set up the very basics */

  /* First get the basic information out of the file. */
  if (read_the_file(fp, manifold) == func_failed)
    return func_failed;

  /* Then add the bells and whistles.  Notice that the
     original peripheral curves are preserved for an
     orientable manifold, but replaced for a nonorientable
     one. */
  orient(manifold);
  orient_edge_classes(manifold);
  if ((manifold)->orientability == nonorientable_manifold)
    peripheral_curves(manifold);
  else
    fix_orientation(manifold);
  find_complete_hyperbolic_structure(manifold);
  if ((manifold)->orientability == nonorientable_manifold)
	/*
	 *  peripheral_curves() will have trashed the original
	 *  curves, even on the orientable cusps.
	 *  Repair the damage as best we can.
	 */
    install_shortest_bases(manifold);
  
  int num_read;
  char keyword[100];
  
  double csu, csp; 
  num_read = fscanf(fp,
		    "%99s %lf %lf",
		    keyword, &csu, &csp); 

  (manifold)->CS_value[ultimate] = csu; 
  (manifold)->CS_value[penultimate] = csp;

  (manifold)->CS_value_is_known = (num_read == 3 && strcmp(keyword, "CS") == 0);

  compute_CS_fudge_from_value(manifold);

  return func_OK;
}

FuncResult read_the_file(
    FILE            *fp,
    Triangulation   *manifold)
{
    int             i,
		    c,
		    cusp_index,
		    vertex,
		    face,
		    edge,
		    count,
		    nbr_index,
		    edge_index,
		    digit,
		    d;
    Tetrahedron     **tal,  /* tetrahedron address list */
		    *tet;
    EdgeClass       **eal,  /* edge class address list  */
		    *ec;
    Cusp            **cal,  /* cusp address list        */
		    *cusp;

    if (fscanf(fp, "%d%d%d%d%*lf%*d",
	    &manifold->num_tetrahedra,
	    &manifold->num_cusps,
	    &manifold->num_nonor_cusps,
	    &manifold->orientability) != 4)
	return(func_failed);

    for (i = 0; i < manifold->num_tetrahedra; i++)
	fscanf(fp, "%*d");  /* ignore edge class sizes */

    manifold->num_or_cusps = manifold->num_cusps - manifold->num_nonor_cusps;

    if (manifold->num_tetrahedra < 1
     || manifold->orientability < 0 || manifold->orientability > 2
     || manifold->num_or_cusps + manifold->num_nonor_cusps != manifold->num_cusps)
	return(func_failed);

    cal = NEW_ARRAY(manifold->num_cusps, Cusp *);
    for (count = 0; count < manifold->num_cusps; count++) {
	cusp = NEW_STRUCT(Cusp);
	initialize_cusp(cusp);
	cusp->index = count;
	INSERT_BEFORE(cusp, &manifold->cusp_list_end);
	cal[count] = cusp;
    }
    for (count = 0; count < manifold->num_nonor_cusps; count++) {
	cusp = cal[count];
	cusp->topology = Klein_cusp;
	cusp->m = 0;
	cusp->l = 0;
	cusp->is_complete = 1;
    }
    for (count = manifold->num_nonor_cusps; count < manifold->num_cusps; count++) {
	cusp = cal[count];
	cusp->topology = torus_cusp;
	cusp->m = 0;
	cusp->l = 0;
	cusp->is_complete = 1;
    }

    tal = NEW_ARRAY(manifold->num_tetrahedra, Tetrahedron *);
    for (count = 0; count < manifold->num_tetrahedra; count++) {
	tet = NEW_STRUCT(Tetrahedron);
	initialize_tetrahedron(tet);
	for (i = 0; i < 2; i++)
	    tet->shape[i] = NEW_STRUCT(TetShape);
	INSERT_BEFORE(tet, &manifold->tet_list_end);
	tal[count] = tet;
    }

    eal = NEW_ARRAY(manifold->num_tetrahedra, EdgeClass *);
    for (count = 0; count < manifold->num_tetrahedra; count++) {
	ec = NEW_STRUCT(EdgeClass);
	initialize_edge_class(ec);
	INSERT_BEFORE(ec, &manifold->edge_list_end);
	eal[count] = ec;
    }

    for (count = 0; count < manifold->num_tetrahedra; count++) {
	tet = tal[count];
	tet->index = count;
	for (face = 0; face < 4; face++)
	    if (fscanf(fp, "%d", &nbr_index) != 1
	     || nbr_index < 0 || nbr_index >= manifold->num_tetrahedra)
		return(func_failed);
	    else
		tet->neighbor[face] = tal[nbr_index];
	for (face = 0; face < 4; face++)
	    for (digit = 4; --digit >= 0; )
		if (fscanf(fp, "%1d", &d) != 1
		 || d < 0 || d > 3)
		    return(func_failed);
		else
		    tet->gluing[face] = (tet->gluing[face] << 2) + d;
	for (vertex = 0; vertex < 4; vertex++) {
	    if (fscanf(fp, "%d", &cusp_index) != 1
	     || cusp_index < 0 || cusp_index >= manifold->num_cusps)
		return(func_failed);
	    tet->cusp[vertex] = cal[cusp_index];
	}
	for (c = 0; c < 2; c++)
	    for (vertex = 0; vertex < 4; vertex++)
		for (face = 0; face < 4; face++) {
		    if (fscanf(fp, "%d", &tet->curve[c][right_handed][vertex][face]) != 1)
			return(func_failed);
		    tet->curve[c][left_handed][vertex][face] = 0;
		}
	for (edge = 0; edge < 6; edge++)
	    if (fscanf(fp, "%d", &edge_index) != 1
	     || edge_index < 0 || edge_index >= manifold->num_tetrahedra)
		return(func_failed);
	    else {
		tet->edge_class[edge] = eal[edge_index];
		tet->edge_class[edge]->order++;
		tet->edge_class[edge]->incident_tet         = tet;
		tet->edge_class[edge]->incident_edge_index  = edge;
	    }
	fscanf(fp, "%*d");  /* ignore tet orientation */
	double re, im; 
	fscanf(fp, "%lf%lf", &re, &im);

	tet->shape[complete]->cwl[ultimate][0].rect.real = re;
	tet->shape[complete]->cwl[ultimate][0].rect.imag = im;

	for (i = 0; i < 2; i++)
	    tet->shape[complete]->cwl[ultimate][i+1].rect
		= complex_div( One, complex_minus(One, tet->shape[complete]->cwl[ultimate][i].rect) );  
	for (i = 0; i < 3; i++)
	    tet->shape[complete]->cwl[ultimate][i].log
		= complex_log(tet->shape[complete]->cwl[ultimate][i].rect, PI/3);
	for (i = 0; i < 3; i++)
	    tet->shape[complete]->cwl[penultimate][i]
		= tet->shape[complete]->cwl[ultimate][i];
	*tet->shape[filled] = *tet->shape[complete];

	tet->extra = NULL;

    }

    my_free_array(tal);
    my_free_array(eal);
    my_free_array(cal);

    return(func_OK);
}


void write_old_manifold_file(
    FILE            *fp,
    Triangulation   *manifold)
{
    Tetrahedron *tet;
    EdgeClass   *edge;
    int         count;
    Permutation p;
    int         i,
		j,
		c;

    if (manifold->orientability != oriented_manifold)
	uAcknowledge("Can't save peripheral curves of nonorientable cusps correctly in the old file format.");

    /*
     *  Put indices on the Tetrahedra and EdgeClasses.
     *  The cusps should already have them.
     */
    count = 0;
    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet->next)

	tet->index = count++;

    count = 0;
    for (edge = manifold->edge_list_begin.next;
	 edge != &manifold->edge_list_end;
	 edge = edge->next)

	edge->index = count++;

    /*
     *  Write the new format manifold in the old file format.
     */

    fprintf(fp, "%d\n%d %d\n%d\n\n%.15lf\n15\n",
	manifold->num_tetrahedra,
	manifold->num_cusps,
	manifold->num_nonor_cusps,
	(double)manifold->orientability,
	(double)volume(manifold, NULL));

    for (edge = manifold->edge_list_begin.next;
	 edge != &manifold->edge_list_end;
	 edge = edge->next)

	fprintf(fp, "%d ", edge->order);

    fprintf(fp, "\n");

    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet->next)
    {
	for (i = 0; i < 4; i++)
	    fprintf(fp, " %3d ", tet->neighbor[i]->index);
	fprintf(fp, "\n");

	for (i = 0; i < 4; i++) {
	    p = tet->gluing[i];
	    fprintf(fp, "%d%d%d%d ", p>>6, (p>>4)&3, (p>>2)&3, p&3);
	}
	fprintf(fp, "\n");

	for (i = 0; i < 4; i++)
	    fprintf(fp, " %3d ", tet->cusp[i]->index);
	fprintf(fp, "\n");

	for (c = 0; c < 2; c++) {
	    for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		    fprintf(fp, " %3d ", tet->curve[c][0][i][j] + tet->curve[c][1][i][j]);
	    fprintf(fp, "\n");
	}

	for (i = 0; i < 6; i++)
	    fprintf(fp, " %3d ", tet->edge_class[i]->index);
	fprintf(fp, "\n");

	fprintf(fp, "0\n");

	fprintf(fp, "%.15lf %.15lf\n", 
	    (double)tet->shape[complete]->cwl[ultimate][0].rect.real, 
	    (double)tet->shape[complete]->cwl[ultimate][0].rect.imag);

	fprintf(fp, "\n");
    }

    if (manifold->CS_value_is_known == TRUE)
	fprintf(fp,
		"CS %.18lf %.18lf\n",
		(double)manifold->CS_value[ultimate],
		(double)manifold->CS_value[penultimate]);

    return;
}



/*
 *  At first I thought I just needed to make sure the
 *  peripheral curves adhered to the orientation convention,
 *  but I now realize that the whole orientation of the
 *  manifold is wrong.
 */

#define FIX_UP_VERSION 1

static void fix_orientation(
    Triangulation *manifold)
{
    copy_curves_to_scratch(manifold, 0, FALSE);
    copy_curves_to_scratch(manifold, 1, FALSE);

    compute_intersection_numbers(manifold);

#if (FIX_UP_VERSION == 0)
    reverse_meridians_where_necessary(manifold);
#endif

#if (FIX_UP_VERSION == 1)
{
    int         standard;
    Cusp        *cusp;
    Tetrahedron *tet;
    int         i, j, k;

    /*
     *  All the intersection numbers ought to be the same.
     */
    standard = manifold->cusp_list_begin.next->intersection_number[L][M];
    for (cusp = manifold->cusp_list_begin.next;
	 cusp != &manifold->cusp_list_end;
	 cusp = cusp->next)
	if (cusp->intersection_number[L][M] != standard)
	    uFatalError("fix_orientation", "compatibility");

    /*
     *  If the curves are backwards, reorient the manifold.
     */
    if (standard == -1) {

	reorient(manifold);

	/*
	 *  Now the manifold's fine, but reorient() has reversed the
	 *  meridians, on the assumption that they were correctly
	 *  oriented to begin with, which they weren't.  Fix 'em.
	 */
	for (tet = manifold->tet_list_begin.next;
	     tet != &manifold->tet_list_end;
	     tet = tet->next)
	    for (i = 0; i < 2; i++)
		for (j = 0; j < 4; j++)
		    for (k = 0; k < 4; k++)
			tet->curve[M][i][j][k] = - tet->curve[M][i][j][k];
    }
}
#endif
}


#if (FIX_UP_VERSION == 0)

static void reverse_meridians_where_necessary(
    Triangulation   *manifold)
{
    Tetrahedron *tet;
    int         i,
		j,
		k;

    /* which Tetrahedron */
    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet->next)

	/* which ideal vertex */
	for (i = 0; i < 4; i++)

	    if (tet->cusp[i]->intersection_number[L][M] == -1)

		/* which side of the vertex */
		for (j = 0; j < 4; j++)

		    if (i != j)

			/* which sheet (right_handed or left_handed) */
			for (k = 0; k < 2; k++)

			    tet->curve[M][k][i][j] = - tet->curve[M][k][i][j];

    return;
}

#endif

/*
 *  rehydrate_census_manifold.c  (unix version written by cdh)
 *
 *  This file, which is intended for use in auxiliary programs
 *  rather than SnapPea itself, rehydrates a census manifold
 *  from the TRS# resources.
 */

Triangulation* get_census_manifold(FILE* fp, int census, long n)
{
  Triangulation* theTriangulation = NEW_STRUCT(Triangulation); 
  TersestTriangulation theTersest;

  fseek(fp, n * 18L, 0); /* find start of data */
  fread(theTersest, 17, 1, fp); /* read 17 bytes */
  initialize_triangulation(theTriangulation); /* set up the very basics */
  rehydrate_census_manifold(theTersest, census, n, theTriangulation);

  return theTriangulation; 
}

double census_manifold_volume(FILE* fp, int census, long n)  
{
  Triangulation* theTriangulation = get_census_manifold(fp, census, n); 
  if (!theTriangulation) return 0.; 

  double vol = volume(theTriangulation, 0); 

  free_triangulation(theTriangulation); 
  return vol; 
}

bool locate_volume(FILE* fp, double vol, double eps, int census, int& start)
{
  if (census < 5 || census > 7) return false; 
  double max_vol = REG_TET_VOL * census; 

  if (vol < 0. || vol > max_vol+eps) return false;

  int lo = -1; 
  int hi = num_census_manifolds[census]; 
  int mid; 

  double mid_vol, hi_vol=1000.; 

  while (lo < hi-1) {

    mid = (lo + hi)/2;
    mid_vol = census_manifold_volume(fp, census, mid); 

    if (mid_vol+eps > vol) { hi=mid; hi_vol=mid_vol; }
    else { lo=mid; }
  }
  
  if (hi_vol-eps < vol) {
    start = hi; 
    return true; 
  }

  return false; 
}

bool is_census_manifold(FILE* fp, Triangulation* T, double eps, int census, int& n)
{
  double vol = volume(T, 0); 
  if (!locate_volume(fp, vol, eps, census, n)) return false; 

  Triangulation* Tn;
  double Tn_vol; 
  int nhi = num_census_manifolds[census];
  Boolean isometric=FALSE; 

  while (n < nhi) {
    Tn = get_census_manifold(fp, census, n); 
    Tn_vol = volume(Tn, 0); 
    if (Tn_vol > vol+eps) break; 
    
    if (compute_isometries(T, Tn, &isometric, NULL, NULL)!=func_OK) {
      printf("Problem checking if manifolds are isometric!\n");
    } else if (isometric) break; 

    free_triangulation(Tn); 
    n++;
  }

  if (n < nhi) free_triangulation(Tn); 

  return (isometric==TRUE); 
}

bool find_census_manifold(const char* path, Triangulation* T, int& census, int& n)
{
  char filename[20];
  strcpy(filename, "trs");
  FILE* fp;

  for (census=5; census<8; census++) {
    filename[3] = '0'+ census;       /*files: trs5, trs6,trs7 */
    filename[4] = 0;

    fp = locate_file(path, filename, "rb"); /* read binary file */
    if (!fp) {
      printf("Unable to find file %s\n", filename); 
      return false; 
    }
    
    if (is_census_manifold(fp, T, 1e-6, census, n)) {
      fclose(fp); 
      return true; 
    }
    fclose(fp); 
  }
  return false;
}


static FuncResult read_census_manifold(
    const char *directory,                                 
    int which_census,
    int which_manifold,
    Triangulation *theTriangulation)
{
    TersestTriangulation    theTersest;
    char                    theName[6];
    FILE                    *fp;
    char                    filename[20];
    char                    path_name[100]; 

    if (which_census < 5 || which_census > 7 || which_manifold < 0 ||
	which_manifold >= num_census_manifolds[which_census]) {
      printf("census manifold %d %d does not exist\n", which_census, which_manifold); 
      return func_failed; 
    }

    switch (which_census) {

	case 5:
	    theName[0] = 'm';
	    theName[1] = '0' + (which_manifold / 100) % 10;
	    theName[2] = '0' + (which_manifold /  10) % 10;
	    theName[3] = '0' + (which_manifold /   1) % 10;
	    theName[4] = 0;
	    break;

	case 6:
	    theName[0] = 's';
	    theName[1] = '0' + (which_manifold / 100) % 10;
	    theName[2] = '0' + (which_manifold /  10) % 10;
	    theName[3] = '0' + (which_manifold /   1) % 10;
	    theName[4] = 0;
	    break;

	case 7:
	    theName[0] = 'v';
	    theName[1] = '0' + (which_manifold / 1000) % 10;
	    theName[2] = '0' + (which_manifold /  100) % 10;
	    theName[3] = '0' + (which_manifold /   10) % 10;
	    theName[4] = '0' + (which_manifold /    1) % 10;
	    theName[5] = 0;
	    break;

	default:
	    return func_failed; 
    }

    strcpy(filename,"trs");
    filename[3] = '0'+ which_census;       /*files: trs5, trs6,trs7 */
    filename[4] = 0;

    fp = locate_file(directory, filename,"rb"); /* read binary file */
    if (!fp) {
      printf("Unable to find file %s\n", filename); 
      return func_failed; 
    }

    fseek(fp,(long) which_manifold*18,0); /* find start of data */
    fread(theTersest, 17, 1, fp); /* read 17 bytes */
    fclose(fp);

    initialize_triangulation(theTriangulation); /* set up the very basics */

    rehydrate_census_manifold(theTersest, which_census, which_manifold, theTriangulation);

    set_triangulation_name(theTriangulation, theName);

    return func_OK;
}

Triangulation* read_census_manifold(const char* path, int which_census, int which_manifold)
{
  Triangulation* manifold = NEW_STRUCT(Triangulation); 
  if (!manifold) return NULL;

  if (read_census_manifold(path, which_census, which_manifold, manifold)==func_OK)
    return manifold;

  my_free(manifold);
  return NULL;
}

/* Quick and dirty implementation of some UI functions
for use in testing the kernel. */

void uAcknowledge(
    char *message)
{
    printf("%s\n", message);
    return;
}

void uFatalError(
    char    *function,
    char    *file)
{
    printf("A fatal error has occurred in the function %s() in the file %s.c.\n",
	function, file);
    exit(0);
}

void uLongComputationBegins(
    char    *message,
    Boolean is_abortable)
{
    printf("%s\n", message);
    return;
}

FuncResult uLongComputationContinues()
{
    return(func_OK);
}

void uLongComputationEnds()
{
    return;
}

/* added by cdh */
/*
 *  Presents the string *message to the user and asks the user to choose
 *  one of the responses.  Returns the number of the chosen response
 *  (numbering starts at 0).  In an interactive context, the UI should
 *  present the possible responses evenhandedly -- none should be
 *  presented as a default.  However, in a batch context (when no human
 *  is present), uQuery should return the default_response.
 */

int uQuery(
    char    *message,
    int     num_responses,
    char    *responses[],
    int     default_response)
{
    int     i,response;
    printf("%s\n",message);
    for(i=0;i<num_responses;i++)
	printf("Type %d to %s\n",i,responses[i]);
    printf("What do you want to do?   ");
    scanf("%d",&response);
    if(0 <= response && response <num_responses) return(response);
	else return(default_response);

}

