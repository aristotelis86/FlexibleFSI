/**********************************************************************
      FlexibleSheet class: Creates the flexible 1D object 
      consisting of control points and springs. Properties
      related to its behaviour are also created.
      Methods to handle and update the object can be found 
      here.
      
Example code:
      (To be filled)

**********************************************************************/

class FlexibleSheet extends LineSegBody {
  //============================= Attributes =================================//
  int numOfpoints;
  int numOfsprings;
  float Length, Mass, segLength, pointMass, stiffness, damping;
  //float [] xcurrent, ycurrent, vxcurrent, vycurrent;
  //PVector [] accelCurrent, accelPred;
  //float [] xPred, yPred, vxPred, vyPred;
  //float [] xnew, ynew, vxnew, vynew;
  float dtmax;
  
  // Dependancy from other objects
  ControlPoint [] cpoints; // stores the control points of the system
  Spring [] springs; // stores the connecting springs
  
  FlexibleSheet( float L_, float th_, float M_, int resol, float stiffness_, PVector lpos, PVector align, Window window ) {
    super(lpos.x, lpos.y, window); 
    thk = th_;
    weight = window.pdx(thk);
    
    align.normalize();
    PVector X0;
    for ( int i=0; i <= L_; i+=resol ) {
      X0 = PVector.add( lpos, PVector.mult(align, i));
      super.add(X0.x, X0.y);
    }
    super.end(false);
    
    Length = this.coords.get(0).dist(this.coords.get(this.coords.size()-1));
    pointMass = M_;
    stiffness = stiffness_;
    numOfpoints = this.coords.size();
    numOfsprings = numOfpoints - 1;
    segLength = resol; // resting length of each spring
    Mass = pointMass * numOfpoints;
    //damping = Determine_damping();
    damping = 0;
    
    cpoints = new ControlPoint[numOfpoints];
    springs = new Spring[numOfsprings];
    
    for (int i = 0; i < numOfpoints; i++) cpoints[i] = new ControlPoint( this.coords.get(i), pointMass, thk, window );
    for (int i = 0; i < numOfsprings; i++) springs[i] = new Spring( cpoints[i], cpoints[i+1], segLength, stiffness, damping, thk, window );
    
    //dtmax = Determine_time_step();
    dtmax = 0.005;
    //InitUpdateVars();
  } // end of Constructor
  
  //FlexibleSheet( float L_, float th_, float M_, float stiffness_, float x, float y, PVector align, Window window ) {
  //  this( L_, th_, M_, 1, stiffness_, x, y, align, window);
  //}
  
  //============================= Methods =================================//
  // !!!! override parent class method !!!!
  // Provides info to the fluid solver that it's not stationary object 
  boolean unsteady() {return true;}
  
  // !!!! override parent class method !!!!
  // Velocity of the fluid at some point from the surface of the object 
  float velocity( int d, float dt, float x, float y ) {
    int coorSize = super.coords.size();
    float sumWeights = 0;
    float sumVel = 0;
    
    for (int j=0; j<coorSize; j++) {
      PVector r = new PVector(x, y);
      r.sub(coords.get(j));
      float factor = exp((-1)*r.magSq());
      
      if (d==1) {
        sumVel += cpoints[j].velocity.x * factor;
        sumWeights += factor;
      }
      else {
        sumVel += cpoints[j].velocity.y * factor;
        sumWeights += factor;
      }
    }
    if (sumWeights<1e-9) return 0.0;
    else return sumVel/sumWeights;
  } // end velocity()
  
  // Display control points for testing
  void pointsDisplay() {
    for (ControlPoint cp : cpoints) cp.display();
  }
  
  // Calculate current length of sheet
  float CurrentLength() {
    float newL = 0;
    for (int i=0; i<numOfpoints-1; i++) newL += this.coords.get(i).dist(this.coords.get(i+1));
    return newL;
  }
  
  //// Initialize state arrays for updating
  //void InitUpdateVars() {
  //  xcurrent = new float[numOfpoints]; xPred = new float[numOfpoints]; xnew = new float[numOfpoints];
  //  ycurrent = new float[numOfpoints]; yPred = new float[numOfpoints]; ynew = new float[numOfpoints];
  //  vxcurrent = new float[numOfpoints]; vxPred = new float[numOfpoints]; vxnew = new float[numOfpoints];
  //  vycurrent = new float[numOfpoints]; vyPred = new float[numOfpoints]; vynew = new float[numOfpoints];
  //  accelCurrent = new PVector[numOfpoints];
  //  accelPred = new PVector[numOfpoints];
  //}
  
  //// Calculate damping coefficient for numerical purposes
  //float Determine_damping() {
  //  float d = 2*sqrt(stiffness*pointMass);
  //  return d;
  //}
  
  //// Determine the size of the time step
  //float Determine_time_step() {
  //  float ReLam, ImLam, dt;
  //  float n = float(numOfsprings);
  //  //float n = float(numOfpoints);
    
  //  float RootDet = 1-stiffness*pointMass*sq(Length)/(2*sq(damping)*sq(n));
  //  float fact = 4*damping*sq(n)/(pointMass*sq(Length));
    
  //  if (RootDet<0) {
  //    ReLam = -fact;
  //    ImLam = -fact*sqrt(-RootDet);
  //  }
  //  else {
  //    ReLam = fact*(-1-sqrt(RootDet));
  //    ImLam = 0.0;
  //  }
    
  //  dt = -1*ReLam/(sq(ReLam)+sq(ImLam));
  //  //dt = -1/(sq(ReLam)+sq(ImLam));
  //  return dt;
  //} // end of Determine_time_step
  
  //// Determine relative positions based on the stiffness
  //void Calculate_Stretched_Positions( PVector F ) {
    
  //  float g = F.mag();
  //  float ll;
  //  PVector newll;
  //  float [] xnew = new float[numOfpoints];
  //  float [] ynew = new float[numOfpoints];
  //  PVector FDir = F.copy();
  //  FDir.normalize();
    
  //  xnew[0] = cpoints[0].position.x; ynew[0] = cpoints[0].position.y;
  //  for (int i = 0; i < numOfsprings; i++) {
  //    ll = (float(numOfpoints) - float(i) - 1)*(pointMass*g/stiffness) + segLength;
  //    newll = PVector.add(cpoints[i].position, PVector.mult(FDir,ll));
  //    cpoints[i+1].UpdatePosition(newll.x, newll.y);
  //  }
  //} // end of Calculate_Stretched_Positions
  
  //// Update the state of the sheet 
  //void UpdateState( float [] Xnew, float [] Ynew, float [] VXnew, float [] VYnew ) {
  //  for (int i = 0; i < numOfpoints; i++) {
  //    cpoints[i].UpdatePosition(Xnew[i], Ynew[i]);
  //    cpoints[i].UpdateVelocity(VXnew[i], VYnew[i]);
  //  }
  //  getOrth();
  //  getBox();
  //}
  
  //// Update the state of the sheet 
  //void UpdateState( float [] Xnew, float [] Ynew ) {
  //  for (int i = 0; i < numOfpoints; i++) {
  //    cpoints[i].UpdatePosition(Xnew[i], Ynew[i]);
  //  }
  //  getOrth();
  //  getBox();
  //}
  
  //// Get current positions and velocities
  //void getState() {
  //  for(int i=0; i<numOfpoints; i++) {
  //    xcurrent[i] = cpoints[i].position.x;
  //    ycurrent[i] = cpoints[i].position.y;
  //    vxcurrent[i] = cpoints[i].velocity.x;
  //    vycurrent[i] = cpoints[i].velocity.y;
  //  }
  //}
  
  //// Move for testing purposes (translation only)
  //void move() {
  //  for (ControlPoint cp : cpoints) {
  //    cp.velocity = new PVector(0.1, .2);
  //    cp.position.add( cp.velocity );
  //  }
  //  getOrth();
  //  getBox();
  //}
  
  //// Prescribed motion of the sheet to test effect on flow
  //void waveOnSheet( float t ) {
    
  //  float sinAmp = Length/5.;
  //  float sinN = 1.;
  
  //  int nn = cpoints.length;
  //  float [] x = new float[nn];
  //  float [] y = new float[nn];
  //  float [] vx = new float[nn];
  //  float [] vy = new float[nn];
    
  //  x[0] = cpoints[0].position.x;
  //  y[0] = cpoints[0].position.y;
  //  vx[0] = 0.;
  //  vy[0] = 0.;
  //  for (int i = 1; i < nn; i++) {
  //    x[i] = (Length/(nn-1)) + x[i-1];
  //    y[i] = (sinAmp * sin(sinN*PI*(x[i]-x[0])/Length))*sin(2*t) + y[0];
  //    vx[i] = 0.;
  //    vy[i] = 2*cos(2*t)*(sinAmp * sin(sinN*PI*(x[i]-x[0])/Length));
  //  }
  //  UpdateState(x, y, vx, vy);
  //  getOrth();
  //  getBox();
  //}  
  
  //// Apply internal Forces to control points
  //void ApplyIntForces() {
  //  for (Spring s : springs) s.applyAllForces();
  //}
  
  //// Apply external forces to control points
  //void ApplyExtForces( PVector [] F) {
  //  for (int i=0; i<numOfpoints; i++) cpoints[i].force.add(F[i]);
  //}
  
  //// Apply gravitational Forces to particles
  //void ApplyGravity( PVector g) {
  //  for (ControlPoint cp : cpoints) {
  //    PVector extg = g.copy();
  //    extg.mult(cp.mass);
  //    cp.applyForce(extg);
  //  }
  //}
  
  //// Clear all forces acting in control points
  //void ClearForces() {
  //  for (ControlPoint p : cpoints) p.clearForce();
  //}
  
  //// calculate the pressure force at each point
  //PVector [] pressForcePoints ( Field p ) {
    
  //  int orthSize = orth.length;
  //  PVector [] pf = new PVector[numOfpoints];
  //  for (int i=0; i<numOfpoints; i++) pf[i] = new PVector(0, 0);
    
  //  for ( int s=-1; s<=1; s+=2 ) {
      
  //    for ( int j=0; j<orthSize; j++ ) {
  //      float pdl = p.linear( cpoints[j].position.x+0.5*s*thk*orth[j].nx, cpoints[j].position.y+0.5*s*thk*orth[j].ny )*orth[j].l;
  //      PVector pTemp = new PVector(s*pdl*orth[j].nx, s*pdl*orth[j].ny);
  //      pf[j].sub(pTemp);
  //      pf[j+1].sub(pTemp);
  //    }
  //  }
  //  for (int j=1; j<numOfpoints-1; j++) pf[j].div(2);
    
  //  return pf;
  //}
  
  //// calculate the pressure force at each point
  //PVector [] stressForcePoints ( VectorField Vel, float nu ) {
    
  //  Field omega = Vel.curl();
    
  //  int orthSize = orth.length;
  //  PVector [] ps = new PVector[numOfpoints];
  //  for (int i=0; i<numOfpoints; i++) ps[i] = new PVector(0, 0);
    
  //  for ( int s=-1; s<=1; s+=2 ) {
  //    for ( int j=0; j<orthSize; j++ ) {
  //      float omegaZVal = omega.linear( cpoints[j].position.x+0.5*s*thk*orth[j].nx, cpoints[j].position.y+0.5*s*thk*orth[j].ny );
  //      PVector pTemp = new PVector(s*omegaZVal*orth[j].ny, s*omegaZVal*orth[j].nx);
  //      ps[j].add(pTemp);
  //      ps[j+1].add(pTemp);
  //    }
  //  }
  //  for (int j=1; j<numOfpoints-1; j++) ps[j].div(2);
  //  for (int j=0; j<numOfpoints; j++) ps[j].mult(nu);
    
  //  return ps;
  //}
  
  //// Update based on Predictor-Corrector Scheme (gravity only)
  //void update(float dt, PVector p) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  // Apply Forces for this step
  //  ClearForces();
  //  ApplyIntForces();
  //  ApplyGravity( p );
    
  //  // calculate acceleration
  //  for (ControlPoint cp : cpoints) {
  //    if (!cp.fixed) cp.update( dt );
  //  }
  //  getOrth();
  //  getBox();
  //} // end of update (prediction)
  
  //void update2(float dt, PVector p) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  // Apply Forces for the correction
  //  ClearForces();
  //  ApplyIntForces();
  //  ApplyGravity( p );
    
  //  // calculate acceleration for the correction
  //  for (ControlPoint cp : cpoints) {
  //    if (!cp.fixed) cp.update2( dt );
  //  }
  //  getOrth();
  //  getBox();
  //} // end of Trapezoid
  
  //// Alternative Update based on Predictor-Corrector Scheme (gravity only)
  //void updateAlt(float dt, PVector p) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  // Apply Forces for this step
  //  ClearForces();
  //  ApplyIntForces();
  //  ApplyGravity( p );
    
  //  // calculate acceleration
  //  for (ControlPoint cp : cpoints) {
  //    if (!cp.fixed) cp.updateAlt( dt );
  //  }
  //  getOrth();
  //  getBox();
  //} // end of update (prediction)
  
  //void updateAlt2(float dt, PVector p) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  // Apply Forces for the correction
  //  ClearForces();
  //  ApplyIntForces();
  //  ApplyGravity( p );
    
  //  // calculate acceleration for the correction
  //  for (ControlPoint cp : cpoints) {
  //    if (!cp.fixed) cp.updateAlt2( dt );
  //  }
  //  getOrth();
  //  getBox();
  //} // end of Trapezoid
  
  
  //// Update based on Predictor-Corrector Scheme (gravity only)
  //void update( float dt, BDIM flow ) { update( dt, flow, new PVector(0,0)); }
  //void update( float dt, BDIM flow, PVector g) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  int N = numOfpoints;
  //  PVector [] pressForce = new PVector[N];
  //  PVector [] stressForce = new PVector[N];
    
  //  pressForce = pressForcePoints ( flow.p );
  //  stressForce = stressForcePoints( flow.u, flow.nu );
    
  //  // Apply Forces for this step
  //  ClearForces();
  //  ApplyIntForces();
  //  ApplyExtForces( pressForce );
  //  ApplyExtForces( stressForce );
  //  ApplyGravity( g );
    
  //  // calculate acceleration
  //  for (ControlPoint cp : cpoints) {
  //    if (!cp.fixed) cp.update( dt );
  //  }
  //  getOrth();
  //  getBox();
  //} // end of update (prediction)
  
  //void update2( float dt, BDIM flow ) { update2( dt, flow, new PVector(0,0)); }
  //void update2(float dt, BDIM flow, PVector g) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  int N = numOfpoints;
  //  PVector [] pressForce = new PVector[N];
  //  PVector [] stressForce = new PVector[N];
    
  //  pressForce = pressForcePoints ( flow.p );
  //  stressForce = stressForcePoints( flow.u, flow.nu );
    
  //  // Apply Forces for the correction
  //  ClearForces();
  //  ApplyIntForces();
  //  ApplyExtForces( pressForce );
  //  ApplyExtForces( stressForce );
  //  ApplyGravity( g );
    
  //  // calculate acceleration for the correction
  //  for (ControlPoint cp : cpoints) {
  //    if (!cp.fixed) cp.update2( dt );
  //  }
  //  getOrth();
  //  getBox();
  //} // end of Trapezoid
  
  //// Alternative Update based on Predictor-Corrector Scheme (gravity only)
  //void updateAlt( float dt, BDIM flow ) { updateAlt( dt, flow, new PVector(0,0)); }
  //void updateAlt( float dt, BDIM flow, PVector g) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  int N = numOfpoints;
  //  PVector [] pressForce = new PVector[N];
  //  PVector [] stressForce = new PVector[N];
    
  //  pressForce = pressForcePoints ( flow.p );
  //  stressForce = stressForcePoints( flow.u, flow.nu );
    
  //  // Apply Forces for this step
  //  ClearForces();
  //  ApplyIntForces();
  //  ApplyExtForces( pressForce );
  //  ApplyExtForces( stressForce );
  //  ApplyGravity( g );
    
  //  // calculate acceleration
  //  for (ControlPoint cp : cpoints) {
  //    if (!cp.fixed) cp.updateAlt( dt );
  //  }
  //  getOrth();
  //  getBox();
  //} // end of update (prediction)
  
  //void updateAlt2( float dt, BDIM flow ) { updateAlt2( dt, flow, new PVector(0,0)); }
  //void updateAlt2(float dt, BDIM flow, PVector g) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  int N = numOfpoints;
  //  PVector [] pressForce = new PVector[N];
  //  PVector [] stressForce = new PVector[N];
    
  //  pressForce = pressForcePoints ( flow.p );
  //  stressForce = stressForcePoints( flow.u, flow.nu );
    
  //  // Apply Forces for the correction
  //  ClearForces();
  //  ApplyIntForces();
  //  ApplyExtForces( pressForce );
  //  ApplyExtForces( stressForce );
  //  ApplyGravity( g );
    
  //  // calculate acceleration for the correction
  //  for (ControlPoint cp : cpoints) {
  //    if (!cp.fixed) cp.updateAlt2( dt );
  //  }
  //  getOrth();
  //  getBox();
  //} // end of Trapezoid
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  //// ****************************************************************** //
  //// Trapezoid (Predictor-Corrector) Scheme
  //void update(float dt, VectorField Vel, float nu ) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  int N = numOfpoints;
  //  PVector [] stressForce = new PVector[N];
  //  float pMass = pointMass;
    
  //  getState();
  //  storeOldState();
  //  stressForce = stressForcePoints( Vel, nu );
  //  println(stressForce[N-1]);
  
  //  // Apply Forces for this step
  //  ApplyAllForces();
    
  //  // calculate acceleration
  //  for (int i = 0; i < N; i++) {
  //    accelCurrent[i] = cpoints[i].force.copy();
  //    accelCurrent[i].div(pMass);
  //  }
  //  // accumulate any acceleration due to external forces
  //  for (int i = 0; i < N; i++) {
      
  //    PVector ext = stressForce[i].copy();
  //    ext.div(pMass);
  //    accelCurrent[i].add(ext);
  //  }
    
  //  // Calculate estimation
  //  for (int i = 0; i < N; i++) {
  //    if (!cpoints[i].fixed) {
  //      xPred[i] = xcurrent[i] + dt*vxcurrent[i];
  //      yPred[i] = ycurrent[i] + dt*vycurrent[i];
  //      vxPred[i] = vxcurrent[i] + dt*accelCurrent[i].x;
  //      vyPred[i] = vycurrent[i] + dt*accelCurrent[i].y;
  //    }
  //    else {
  //      xPred[i] = xcurrent[i];
  //      yPred[i] = ycurrent[i];
  //      vxPred[i] = vxcurrent[i];
  //      vyPred[i] = vycurrent[i];
  //    }
  //  }
  //  // Update the state of the filament for the correction
  //  UpdateState(xPred, yPred, vxPred, vyPred);
  //} // end of update (prediction)
  
  
  //void update2(float dt, VectorField Vel, float nu ) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  int N = numOfpoints;
  //  PVector [] stressForce = new PVector[N];
  //  float pMass = pointMass;
    
  //  stressForce = stressForcePoints( Vel, nu );
    
  //  // Apply Forces for the correction
  //  ApplyAllForces();
    
  //  // calculate acceleration for the correction
  //  for (int i = 0; i < N; i++) {
  //    accelPred[i] = cpoints[i].force.copy();
  //    accelPred[i].div(pMass);
  //  }
  //  // accumulate any acceleration due to external forces
  //  for (int i = 0; i < N; i++) {
  //    PVector ext = stressForce[i].copy();
  //    ext.div(pMass);
  //    accelPred[i].add(ext);
  //  }
    
  //  // Calculate at the new state
  //  for (int i = 0; i < N; i++) {
  //    if (!cpoints[i].fixed) {
  //      xnew[i] = xcurrent[i] + 0.5*dt*(vxcurrent[i]+vxPred[i]);
  //      ynew[i] = ycurrent[i] + 0.5*dt*(vycurrent[i]+vyPred[i]);
  //      vxnew[i] = vxcurrent[i] + 0.5*dt*(accelCurrent[i].x+accelPred[i].x);
  //      vynew[i] = vycurrent[i] + 0.5*dt*(accelCurrent[i].y+accelPred[i].y);
  //    }
  //    else {
  //      xnew[i] = xcurrent[i];
  //      ynew[i] = ycurrent[i];
  //      vxnew[i] = vxcurrent[i];
  //      vynew[i] = vycurrent[i];
  //    }
  //  }
  //  // Update the state of the filament for the correction
  //  UpdateState(xnew, ynew, vxnew, vynew);
  //} // end of Trapezoid
  
  //// *************************************************************** //
  
  //// Trapezoid (Predictor-Corrector) Scheme
  //void update(float dt, Field p) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  int N = numOfpoints;
  //  PVector [] pressForce = new PVector[N];
  //  float pMass = pointMass;
    
  //  getState();
  //  storeOldState();
  //  pressForce = pressForcePoints ( p );
    
  //  // Apply Forces for this step
  //  ApplyAllForces();
    
  //  // calculate acceleration
  //  for (int i = 0; i < N; i++) {
  //    accelCurrent[i] = cpoints[i].force.copy();
  //    accelCurrent[i].div(pMass);
  //  }
  //  // accumulate any acceleration due to external forces
  //  for (int i = 0; i < N; i++) {
  //    PVector ext1 = pressForce[i].copy();
  //    ext1.div(pMass);
  //    accelCurrent[i].add(ext1);
  //  }
    
  //  // Calculate estimation
  //  for (int i = 0; i < N; i++) {
  //    if (!cpoints[i].fixed) {
  //      xPred[i] = xcurrent[i] + dt*vxcurrent[i];
  //      yPred[i] = ycurrent[i] + dt*vycurrent[i];
  //      vxPred[i] = vxcurrent[i] + dt*accelCurrent[i].x;
  //      vyPred[i] = vycurrent[i] + dt*accelCurrent[i].y;
  //    }
  //    else {
  //      xPred[i] = xcurrent[i];
  //      yPred[i] = ycurrent[i];
  //      vxPred[i] = vxcurrent[i];
  //      vyPred[i] = vycurrent[i];
  //    }
  //  }
  //  // Update the state of the filament for the correction
  //  UpdateState(xPred, yPred, vxPred, vyPred);
  //} // end of update (prediction)

  
  //void update2(float dt, Field p) {
    
  //  if (dt>dtmax) {
  //    println("WARNING dt constrained to maximum permitted:"+dtmax);
  //    dt = dtmax;
  //  }
    
  //  int N = numOfpoints;
  //  PVector [] pressForce = new PVector[N];
  //  float pMass = pointMass;
    
  //  pressForce = pressForcePoints ( p );
    
  //  // Apply Forces for the correction
  //  ApplyAllForces();
    
  //  // calculate acceleration for the correction
  //  for (int i = 0; i < N; i++) {
  //    accelPred[i] = cpoints[i].force.copy();
  //    accelPred[i].div(pMass);
  //  }
  //  // accumulate any acceleration due to external forces
  //  for (int i = 0; i < N; i++) {
  //    PVector ext1 = pressForce[i].copy();
  //    ext1.div(pMass);
  //    accelPred[i].add(ext1);
  //  }
    
  //  // Calculate at the new state
  //  for (int i = 0; i < N; i++) {
  //    if (!cpoints[i].fixed) {
  //      xnew[i] = xcurrent[i] + 0.5*dt*(vxcurrent[i]+vxPred[i]);
  //      ynew[i] = ycurrent[i] + 0.5*dt*(vycurrent[i]+vyPred[i]);
  //      vxnew[i] = vxcurrent[i] + 0.5*dt*(accelCurrent[i].x+accelPred[i].x);
  //      vynew[i] = vycurrent[i] + 0.5*dt*(accelCurrent[i].y+accelPred[i].y);
  //    }
  //    else {
  //      xnew[i] = xcurrent[i];
  //      ynew[i] = ycurrent[i];
  //      vxnew[i] = vxcurrent[i];
  //      vynew[i] = vycurrent[i];
  //    }
  //  }
  //  // Update the state of the filament for the correction
  //  UpdateState(xnew, ynew, vxnew, vynew);
  //} // end of Trapezoid
  
  // *************************************************************** //
  
} //=========== end of FlexibleSheet class ===============