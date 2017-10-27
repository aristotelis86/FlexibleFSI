/**********************************************************************
      ControlPoint class: Creates the control points of 
      a flexible structure and handles their behaviour
      during the runs. 

Example code:
int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution
int N = 3; // number of control points to create 

Window view; // convert pixels to non-dim frame
ControlPoint [] cpoints = new ControlPoint[N];

void settings(){
    size(600, 600);
}

void setup() {
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  for (int i=0; i<N; i++) {
    cpoints[i] = new ControlPoint( new PVector(random(view.x.inE),random(view.y.inE)), 5,  5, view );
  }
}

void draw() {
  
  for (ControlPoint cp : cpoints) {
    cp.randomWalk();
    cp.display();
  }
}

**********************************************************************/

class ControlPoint {
  //================= Attributes ====================//
  // Translation
  PVector position; // current position
  PVector positionOld;
  PVector velocity; // current velocity
  PVector velocityOld;
  PVector acceleration; // current acceleration
  PVector accelerationOld; 
  boolean xfixed; // fix the particle at its y-axis
  boolean yfixed; // fix the particle at its x-ayis
  
  // Rotation
  float theta; // heading
  float omega; // rotational velocity
  float rotAccel; // rotational acceleration
  boolean rfixed; // constrain rotaion
  
  // Dynamics
  PVector force; // force acting on the point-mass
  PVector impForce; // impact force tracking
  float mass; // mass of the point
  
  // Display
  Window myWindow; // viewing window
  color c; // for displaying
  float radius; // for collisions model
  
  
  
  //float diameter = 0; // must be removed!!!
  
  
  
  //================= Constructor ====================//
  ControlPoint(PVector position_, float m,  float rad, Window myWindow_) {
    // Translation
    position = position_;
    velocity = new PVector(0, 0);
    acceleration = new PVector(0, 0);
    // Rotation 
    theta = 0;
    omega = 0;
    rotAccel = 0;
    
    // Dynamics
    force = new PVector(0, 0);
    mass = m;
    
    xfixed = false;
    yfixed = false;
    rfixed = false;
    
    // Display
    myWindow = myWindow_;
    c = color(random(1,255), random(1,255), random(1,255));
    radius = rad;
    
    positionOld = position.copy();
    velocityOld = velocity.copy();
    accelerationOld = acceleration.copy();
  }
  
  ControlPoint( PVector position_, float m, Window myWindow_ ) {this(position_, m,  1, myWindow_);}
  
  //================= Methods =====================//
  
  // Display
  void display() {
    noStroke();
    fill(c);
    ellipse(myWindow.px(position.x), myWindow.py(position.y), 2*radius*myWindow.x.r, 2*radius*myWindow.y.r);
    
    pushMatrix();
    translate(myWindow.px(position.x), myWindow.py(position.y));
    rotate(theta);
    stroke(100);
    line(0, 0, 1.2*radius*myWindow.x.r, 0);
    popMatrix();
  }
  
  // Display impact 
  void impDisplay() {
    noStroke();
    fill(255, 0, 0);
    ellipse(myWindow.px(position.x), myWindow.py(position.y), 0.8*2*radius*myWindow.x.r, 0.8*2*radius*myWindow.y.r);
  }
  
  // Update the position
  void setPosition(float x, float y) {
    position.x = x;
    position.y = y;
  }
  
  // Update the velocity
  void setVelocity(float x, float y) {
    velocity.x = x;
    velocity.y = y;
  }
  
  // Store old state
  void StoreOld() {
    positionOld = position.copy();
    velocityOld = velocity.copy();
    accelerationOld = acceleration.copy();
  }
  
  // For testing...
  void randomWalk() {
    this.StoreOld();
    float ex, ey, vx, vy, x, y;
    ex = random(1);
    ey = random(1);
    
    if (ex<.05) vx = -1*velocity.x;
    else vx = velocity.x;
    if (ey<.05) vy = -1*velocity.y;
    else vy = velocity.y;
    
    this.setVelocity(vx,vy);
    x = positionOld.x + velocity.x;
    y = positionOld.y + velocity.y;
    this.setPosition(x,y);
  }
  
  // Clear any forces acting on the particle
  void clearForce() { force.mult(0); }
  
  // Accumulate all the forces acting on the particle
  void ApplyForce(PVector FF) { force.add( FF ); }
  
  // Find the acceleration due to forces
  void calculateAcceleration() {
    PVector accel = force.copy();
    acceleration = accel.div(mass);
  }
  
  // Make the particle free of constraints
  void makeFree() {
    xfixed = false;
    yfixed = false;
    rfixed = false;
  }
  void makeFreex() { xfixed = false; }
  void makeFreey() { yfixed = false; }
  void makeFreer() { rfixed = false; }
  
  // Constrain the particle at its location
  void makeFixed() {
    xfixed = true;
    yfixed = true;
    rfixed = true;
  }
  void makeFixedx() { xfixed = true; }
  void makeFixedy() { yfixed = true; }
  void makeFixedr() { rfixed = true; }
  
  boolean isFixedx() { return xfixed; }
  boolean isFixedy() { return yfixed; }
  boolean isFixedr() { return rfixed; }
  
  // Update methods based on Predictor-Corrector scheme 
  void update( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      StoreOld();
      float x, y, vx, vy;
      if (!this.isFixedx()) {
        x = position.x + t*velocity.x;
        vx = velocity.x + t*acceleration.x;
      }
      else {
        x = position.x;
        vx = velocity.x;
      }
      if (!this.isFixedy()) {
        y = position.y + t*velocity.y;
        vy = velocity.y + t*acceleration.y;
      }
      else {
        y = position.y;
        vy = velocity.y;
      }
      setPosition( x, y );
      setVelocity( vx, vy );
    }
  }
  
  void update2( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      float x, y, vx, vy;
      if (!this.isFixedx()) {
        x = positionOld.x + .5*t*(velocityOld.x + velocity.x);
        vx = velocityOld.x + .5*t*(accelerationOld.x + acceleration.x);
      }
      else {
        x = position.x;
        vx = velocity.x;
      }
      if (!this.isFixedy()) {
        y = positionOld.y + .5*t*(velocityOld.y + velocity.y);
        vy = velocityOld.y + .5*t*(accelerationOld.y + acceleration.y);
      }
      else {
        y = position.y;
        vy = velocity.y;
      }
      setPosition( x, y );
      setVelocity( vx, vy );
    }
  }
  
  // Alternative update methods based on Predictor-Corrector scheme 
  void updateAlt( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      StoreOld();
      float x, y, vx, vy;
      if (!this.isFixedx()) {
        x = position.x + t*velocity.x + 0.5*acceleration.x*sq(t);
        vx = velocity.x + t*acceleration.x;
      }
      else {
        x = position.x;
        vx = velocity.x;
      }
      if (!this.isFixedy()) {
        y = position.y + t*velocity.y + 0.5*acceleration.y*sq(t);
        vy = velocity.y + t*acceleration.y;
      }
      else {
        y = position.y;
        vy = velocity.y;
      }
      setPosition( x, y );
      setVelocity( vx, vy );
    }
  }
  
  void updateAlt2( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      float vx, vy;
      if (!this.isFixedx()) vx = velocityOld.x + .5*t*(accelerationOld.x + acceleration.x);
      else vx = velocity.x;
      if (!this.isFixedy()) vy = velocityOld.y + .5*t*(accelerationOld.y + acceleration.y);
      else vy = velocity.y;
      setVelocity( vx, vy );
    }
  }
  
  // Export info on the motion of the control points
  void dampInfo( PrintWriter outFile ) { this.dampInfo( outFile, 0 ); }
  void dampInfo( PrintWriter outFile, float t ) {
    outFile.println(t+","+position.x+","+position.y+","+velocity.x+","+velocity.y+","+force.x+","+force.y);
  }
  
  // Get the distance between control points
  float distance(ControlPoint other) {
    float d = this.position.dist(other.position);
    return d;
  }
  
  // Rewind the update
  void rewindPosition( float dt ) {
    float x, y;
    x = position.x - dt*(position.x - positionOld.x);
    y = position.y - dt*(position.y - positionOld.y);
    setPosition(x,y);
  }
  
  // Apply impact forces to main force variable
  void applyImpForce() { force.add( impForce ); }
  
  // Add impact forces
  void addImpForce( PVector f ) { impForce.add(f); }
  
  // Clear impact force
  void clearImpForce() { impForce.mult(0); }
  
} // end of ControlPoint class