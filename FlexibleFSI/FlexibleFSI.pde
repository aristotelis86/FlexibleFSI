//////////////////////// INPUTS Section //////////////////////////////

int nx = (int)pow(2,6); // x-dir resolution
int ny = (int)pow(2,6); // y-dir resolution
float xpos = nx/2.; // x-location of leading point
float ypos = ny/4.; // y-location of leading point

float L = ny/4; // length of filament in grid units
int res = 1; // resolution
float thick = 1; // thickness of filament
float M = 1; // line mass density 

float stiff = 500; // stiffness of each spring used (non-dim)

float t = 0; // time keeping variable
float dt; // time step size

PVector align = new PVector(0,1); // initial alignment of filament
PVector gravity = new PVector(0, 5); // constant gravity

//////////////////////// END of INPUTS Section //////////////////////////////
/////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
////////////////////////// Setup Section ///////////////////////////////////
FlexibleSheet [] filament;
Window view; // convert pixels to non-dim frame
WriteInfo myWriter;
//FlexiCollision collider;
CollisionSolver collider;


void settings(){
    size(800, 600);
}

void setup() {
  
  Window view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  filament = new FlexibleSheet[2];
  
  filament[0] = new FlexibleSheet( L, thick, M, res, stiff, xpos, ypos, align, view);
  filament[0].cpoints[0].makeFixed();
  filament[0].Calculate_Stretched_Positions( gravity );
  
  filament[1] = new FlexibleSheet( nx/2., thick, 5*M, res, 3*stiff, nx/4., 7*ny/10., new PVector(1,0), view);
  filament[1].cpoints[0].makeFixed();
  filament[1].cpoints[filament[1].numOfpoints-1].makeFixed();
  
  float dt0 = filament[0].dtmax;
  float dt1 = filament[1].dtmax;
  
  if (dt0<dt1) {
    dt = dt0;
  }
  else {
    dt = dt1;
  }

  myWriter = new WriteInfo(filament);
  //collider = new FlexiCollision(filament);
  collider = new CollisionSolver(filament, view);
} // end of setup


////////////////////////////////////////////////////////////////////////////
////////////////////////// Draw Section ///////////////////////////////////
void draw() {
  background(185);
  fill(0, 0, 0);
  textSize(32);
  text(t, 10, 30);
  
  // Update
  for (FlexibleSheet fs : filament) {
    fs.update(dt, gravity);
    fs.update2(dt, gravity);
  }
  
  //collider.HandleCollisions();
  collider.SolveCollisions();
  // Display
  for (FlexibleSheet fs : filament) fs.display();
  
  // Write Information to files
  for (int i=0; i<2; i++) myWriter.saveInfoSheet(t, gravity, ny, i);
  
  t += dt;
  if (t>3) filament[0].cpoints[0].makeFree();
  saveFrame("./movie/frame_######.png");
  //if (t>120) terminateRun();
  //noLoop();
}


// Gracefully terminate writing...
void keyPressed() {
  myWriter.closeInfos();
  exit(); // Stops the program 
}
void terminateRun() {
  myWriter.closeInfos();
  exit(); // Stops the program 
}







//int nx = (int)pow(2,6); // x-dir resolution
//int ny = (int)pow(2,6); // y-dir resolution
//float t = 0; // time keeping variable
//float dt = 0.01; // time step size
//PVector gravity = new PVector(10, 10); // constant gravity

//int N = 3;

//Window view; // convert pixels to non-dim frame
//ControlPoint [] cpoints = new ControlPoint[N];
//Spring spring;

//void settings(){
//    size(800, 600);
//}

//void setup() {
  
//  Window view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  
//  cpoints[0] = new ControlPoint( new PVector(nx/2.,ny/2.), 3, 3, view);
//  cpoints[1] = new ControlPoint( new PVector(nx/3.,3*ny/4.), 3, 3, view);
//  cpoints[2] = new ControlPoint( new PVector(3,9*ny/10.), 3, 3, view);
  
//  spring = new Spring( cpoints[0], cpoints[1], 10, 10, 1, 1, view );

//}



//void draw() {
//  background(185);
//  fill(0, 0, 0);
//  textSize(32);
//  text(t, 10, 30);
  
//  for (ControlPoint cp : cpoints) cp.clearForce();
//  spring.applyAllForces();
//  for (ControlPoint cp : cpoints) {
//    cp.applyForce(gravity);
//    cp.update(dt);
//  }
//  for (ControlPoint cp : cpoints) cp.BoundCollision();
//  //for (int i=0; i<N-1; i++) {
//  //  for (int j=i+1; j<N; j++) {
//  //    cpoints[i].FastCPointCPointCollision(cpoints[j]);
//  //  }
//  //}
//  //cpoints[2].CPointSpringCollision( spring );
//  cpoints[2].LineSweepsPoint( spring );
//  spring.display();
//  for (ControlPoint cp : cpoints) cp.display();
//  t += dt;
  
//}





//int nx = (int)pow(2,6); // x-dir resolution
//int ny = (int)pow(2,6); // y-dir resolution
//float t = 0; // time keeping variable
//float dt = 0.01; // time step size
//PVector gravity = new PVector(10, 10); // constant gravity

//int N = 3;

//Window view; // convert pixels to non-dim frame
//ControlPoint [] cpoints = new ControlPoint[N];

//void settings(){
//    size(800, 600);
//}

//void setup() {
  
//  Window view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  
//  cpoints[0] = new ControlPoint( new PVector(nx/2.,ny/2.), 3, 5, view);
//  cpoints[1] = new ControlPoint( new PVector(nx/3.,3*ny/4.), 3, 5, view);
//  cpoints[2] = new ControlPoint( new PVector(3,9*ny/10.), 3, 5, view);

//}



//void draw() {
//  background(185);
//  fill(0, 0, 0);
//  textSize(32);
//  text(t, 10, 30);
  
  
//  for (ControlPoint cp : cpoints) {
//    cp.clearForce();
//    cp.applyForce(gravity);
//    cp.update(dt);
//  }
//  for (ControlPoint cp : cpoints) cp.BoundCollision();
//  for (int i=0; i<N-1; i++) {
//    for (int j=i+1; j<N; j++) {
//      cpoints[i].FastCPointCPointCollision(cpoints[j]);
//    }
//  }
//  for (ControlPoint cp : cpoints) cp.display();
//  t += dt;
  
//}













////////////////////////// INPUTS Section //////////////////////////////

//int nx = (int)pow(2,6); // x-dir resolution
//int ny = (int)pow(2,6); // y-dir resolution
//float xpos = nx/2.; // x-location of leading point
//float ypos = 2.; // y-location of leading point

//float L = ny/4; // length of filament in grid units
//int res = 1; // resolution
//float thick = 1; // thickness of filament
//float M = 1; // line mass density 

//float stiff = 500; // stiffness of each spring used (non-dim)

//float t = 0; // time keeping variable
//float dt; // time step size

//PVector align = new PVector(0,1); // initial alignment of filament
//PVector gravity = new PVector(0, 10); // constant gravity

//float maxVel = 300;

////////////////////////// END of INPUTS Section //////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// Setup Section ///////////////////////////////////
//FlexibleSheet filament;
//Window view; // convert pixels to non-dim frame
//WriteInfo myWriter;
//FlexiCollision collider;


//void settings(){
//    size(800, 600);
//}

//void setup() {
  
//  Window view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  
//  filament = new FlexibleSheet( L, thick, M, res, stiff, xpos, ypos, align, view);
//  filament.cpoints[0].makeFixed();
  
//  filament.Calculate_Stretched_Positions( gravity );
  
//  // Create the distortion
//  int N = filament.numOfpoints;
  
//  // Add an impulse (x-dir) to the particles
//  for (int i = 1; i < N; i++) {
//    filament.cpoints[i].velocity.x += ((i-1)/(N-2)) * maxVel;
//  }

//  dt = filament.dtmax;
  
//  myWriter = new WriteInfo(filament);
//  collider = new FlexiCollision(filament);
//} // end of setup


//////////////////////////////////////////////////////////////////////////////
//////////////////////////// Draw Section ///////////////////////////////////
//void draw() {
//  background(185);
//  fill(0, 0, 0);
//  textSize(32);
//  text(t, 10, 30);
  
//  // Update
//  filament.update(dt, gravity);
//  filament.update2(dt, gravity);
  
//  collider.HandleCollisions();
//  // Display
//  filament.display();
  
//  // Write Information to files
//  myWriter.saveInfoSheet(t, gravity, ny, 0);
  
//  t += dt;
  
//  if (t>120) terminateRun();
//  //noLoop();
//}


//// Gracefully terminate writing...
//void keyPressed() {
//  myWriter.closeInfos();
//  exit(); // Stops the program 
//}
//void terminateRun() {
//  myWriter.closeInfos();
//  exit(); // Stops the program 
//}