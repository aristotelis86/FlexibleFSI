/************************ Example Code *******************/
//int nx = 150; // x-dir resolution
//int ny = 150; // y-dir resolution
//int N = 15;
//PVector gravity = new PVector(6,6);
//Window view; // convert pixels to non-dim frame
//SimpleCollisionTest Demo;

//void settings(){
//    size(600, 600);
//}

//void setup() {
//  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
//  Demo = new SimpleCollisionTest( N, gravity, view );
//} // end of setup

//void draw() {
//  background(185); 
//  Demo.runDemo();
//}

//void keyPressed() {
//  Demo.terminateDemo();
//}
/******************************************************/

class SimpleCollisionTest {
  
  //================= Attributes ====================//
  int N;
  Window view; 
  ControlPoint [] cpoints;
  PrintWriter [] myInfo;
  CollisionHandler collider;
  PVector gravity;
  
  //================= Constructor ====================//
  SimpleCollisionTest( int N_, PVector f, Window view_ ) {
    N = N_;
    view = view_;
    gravity = f;
    
    cpoints = new ControlPoint[N];
    myInfo = new PrintWriter[N];
  
    for (int i=0; i<N; i++) {
      cpoints[i] = new ControlPoint( new PVector(random(view.x.inE),random(view.y.inE)), 5,  10, view );
      myInfo[i] = createWriter("./info/cpoints"+i+".txt");
    }
    collider = new CollisionHandler( cpoints );
  }
  
  //================= Methods =====================//
  
  void runDemo() {
    for (ControlPoint cp : cpoints) {
      cp.clearForce();
      cp.ApplyForce( gravity );
      cp.updateAlt( 0.1 );
      cp.updateAlt2( 0.1 );
    }
    
    collider.HandleCollisions();
    for (ControlPoint cp : cpoints) cp.display();
    for (int i=0; i<N; i++) cpoints[i].dampInfo(myInfo[i]);
  }
  
  void terminateDemo() {
    for (int i=0; i<N; i++) {
      myInfo[i].flush();  // Writes the remaining data to the file
      myInfo[i].close();  // Finishes the file
    }
    exit();  // Stops the program
  }
  
  
}