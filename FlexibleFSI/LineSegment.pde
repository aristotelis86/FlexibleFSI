class LineSegment {
  PVector Start;
  PVector End;
  OrthoNormal orth;
  
  LineSegment( PVector s, PVector e ) {
    Start = s;
    End = e;
    orth = new OrthoNormal( Start, End );
  }
  
}