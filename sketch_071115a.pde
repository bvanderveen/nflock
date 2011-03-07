import fullscreen.*;

class NFlock
{
    PApplet pApplet;
    
    class Bird
    {
        public double[] p, v, p0;
        public Bird(double[] p, double[] v) { this.p = p; this.v = v; }
    }
    
    public int 
        count = 200,
        dimension = 6,
        world_width,
        world_height,
        world_depth = 1280;

    double
        phase = 0,
        phaseRate = .0002,
        visibility = .01,
        cohesion = .001,
        alignment = .001,
        separation = .0003,
        drive = .003,
        current_drive = 0,
        attractor_strength = .002,
        momentum = .01,
        randomness = .002;

    public double[][] attractors;

    double[] 
        init_pmax = new double[] { 1, 1, 1, 1, 1, 1 },
        init_pmin = new double[] { 0, 0, 0, 0, 0, 0 },
        init_vmax = new double[] { .05, .05, .05, 0.1, 0.1, 0.1 }, 
        init_vmin = new double[] { -.05, -.05, -.05, -.1, -.1, -.1 },
        random_min = new double[] { -.02, -.02, -.02, -.01, -.01, -.01},
        random_max = new double[] { .02, .02, .02, .01, .01, .01 };
        

    Bird[] birds;
    
    NFlock initialize(int w, int h)
    {
        world_width = w;
        world_height = h;
        size(world_width, world_height);
        background(#000000);
        
        // fix values
        visibility = Math.pow(visibility,2);
        
        birds = new Bird[count];

        for (int i = 0; i < birds.length; i++)
            birds[i] = new Bird(random(init_pmax, init_pmin),  random(init_vmax, init_vmin));
         
        return this;
    }
    
    void iterate()
    {
      for (int j = 0; j < birds.length; j++)
      {
          solve_bird(birds[j],j);
          draw_bird(birds[j]);
      }
      //println();
      
      fill(color(0,0,0,1));
      rect(0,0,world_width,world_height);
    }
    
    void solve_bird(Bird bird, int index)
    {
        //println("solving");
        double[]
            cohesion_point = new double[dimension],
            alignment_vector = new double[dimension],
            separation_vector = new double[dimension];

        System.arraycopy(bird.p, 0, cohesion_point, 0, dimension);
        System.arraycopy(bird.v, 0, alignment_vector, 0, dimension);

        for (int i = 0; i < count; i++)
        {
            if(i == index)
                continue;
                
            if (distance_squared(bird.p, birds[i].p) > visibility)
                continue;

            //println ("visible");
            for (int j = 0; j < dimension; j++)
                cohesion_point[j] += birds[i].p[j];
            
            // simple sum, normalize later
            for (int j = 0; j < dimension; j++)
                alignment_vector[j] += birds[i].v[j];
            
            for (int j = 0; j < dimension; j++)
                separation_vector[j] += bird.p[j] - birds[i].p[j];
        }
        
        double[] driver = new double[dimension];
        
        
        if (current_drive < .000001)
        {
          driver[0] = .5;
          driver[1] = .5;
        }
        else
        {
          driver[0] = (double)mouseX / (double)world_width;
          driver[1] = (double)mouseY / (double)world_height;
        }
        driver[2] = .5;
        driver[3] = .5;// * Math.sin(phase) + .5;
        driver[4] = .5;// * Math.sin(phase * phase * .8 + 1) + .5;
        driver[5] = .5;// * Math.sin(phase + 2) + .5;
        
        
        double[][] components = new double[7][];
        
        components[0] = scale(cohesion, normalize(difference(cohesion_point, bird.p)));
        components[1] = scale(alignment, normalize(alignment_vector));
        components[2] = scale(separation, normalize(separation_vector));
        components[3] = scale(momentum, normalize(bird.v));
        components[4] = scale(randomness, random(random_min, random_max));
        components[5] = scale(current_drive, normalize(difference(driver, bird.p)));
        components[6] = new double[dimension];
        
        double[] attractor = attractors[index % attractors.length];
        
        //if (distance_squared(bird.p, attractor) < visibility) 
        components[6] = scale(attractor_strength, normalize(difference(attractor, bird.p)));
        
        bird.v = sum(components);
    }

    void draw_bird(Bird bird)
    {
        bird.p0 = bird.p;      
        bird.p = sum(new double[][] { bird.p, bird.v });
       
        
        color c = color((int)bird.p[3] * 255,
          (int)(bird.p[4] * 255),
          (int)(bird.p[5] * 255),
          (int)(bird.p[2] * 255));
        stroke(c);
        
        int p0x = (int)(bird.p0[0] * world_width);
        int px = (int)(bird.p[0] * world_width);
        int p0y = (int)(bird.p0[1] * world_height);
        int py = (int)(bird.p[1] * world_height);
    
        //println("x = " + p0x + " y = " + p0y);
        line(
          (int)p0x < 0 ? 0 : p0x,
          (int)p0y < 0 ? 0 : p0y,
          (int)px < 0 ? 0 : px,
          (int)py < 0 ? 0 : py);
    }

    double[] random(double[] upper, double[] lower)
    {
        double[] r = new double[upper.length];
        for (int i = 0; i < r.length; i++)
            r[i] = lower[i] + Math.random() * (upper[i] - lower[i]);
        return r;
    }

    double distance_squared(double[] a, double[] b)
    {
        double r = 0, t = 0;
        for (int i = 0; i < a.length; i++)
        {
            t = (a[i] - b[i]);
            r += t * t;
        }
        return r;
    }

    double magnitude_squared(double[] a)
    {
        return distance_squared(a, new double[a.length]);
    }

    double[] sum(double[][] a)
    {
        double[] r = new double[a[0].length];
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[0].length; j++)
                r[j] += a[i][j];
        return r;
    }

    double[] difference(double[] a, double[] b)
    {
        double[] r = new double[a.length];
        for (int i = 0; i < a.length; i++)
            r[i] = a[i] - b[i];
        return r;
    }

    double[] scale(double a, double[] x)
    {
        double[] r = new double[x.length];
        for (int i = 0; i < x.length; i++)
            r[i] = a * x[i];
        return r;
    }

    double[] normalize(double[] a)
    {
        // slow. better way to find invsqrt?
        
        double m = magnitude_squared(a);
        
        if(Math.abs(m) < .0000000000001)
          m = .0000000000001;
        return scale(1 / Math.sqrt(m), a);
    }   
}

class PathMapper
{
  public ArrayList<Path> paths;

  public PathMapper()
  { 
    paths = new ArrayList<Path>();
    paths.add(new Path());
  }
  
  public void addPoint(double x, double y) 
  {
    double[] a = new double[2];
    a[0] = x;
    a[1] = y;
    
    paths.get(paths.size() - 1).add(a);
  }
  
  public void addPath()
  {
    paths.get(paths.size() - 1).close();
    paths.add(new Path());
  }
  
  public void done()
  {
    paths.get(paths.size() - 1).close();
  }
  
  class Path
  {
    double len;
    double[] lens;
    ArrayList<double[]> vertices;
    
    public Path()
    {
      vertices = new ArrayList<double[]>();
      lens = new double[10];
    }
    
    public void add(double[] v)
    {
      vertices.add(v);
      
      if (vertices.size() > 1)
        addLen(false);
    }
    
    void addLen(boolean last)
    {
      double[] x0 = vertices.get(vertices.size() - 1);
      double[] x1 = vertices.get(last ? 0 : vertices.size() - 2);
      double d0 = (x1[0] - x0[0]);
      double d1 = (x1[1] - x0[1]);
      int il = vertices.size() - (last ? 1 : 2);
      lens[il] = Math.sqrt(d0 * d0 + d1 * d1);
      len += lens[il];
      
      println("added segment of length " + lens[il]);
    }
    
    public void close()
    {
      addLen(true);
      
      println("closed path. segs = " + vertices.size() + ", len = " + len);
      for(int i = 0; i < lens.length; i++)
      {
        println("lens[" + i + "] = " + lens[i]);
      }
    }
    
    public double[] getPosition(double t)
    {
      double p = t * len;
      double cp = 0;
      int count = vertices.size();
      for (int i = 0; i < count; i++)
      {
        double l = lens[i];
        double test = cp + l;
        if (p < test || i == count - 1)
        {
          double ts = (p - cp) / l;
          //println("i = " + i + ", p = " + p + ", cp = " + cp + ", l = " + l + ", t = " + t + " seg = " + i + ", ts = " + ts);
          return positionOnSegment(i, ts);
        }
        cp += l;
      }
      
      throw new Error("dang.");
    }
    
    double[] positionOnSegment(int i, double t)
    {
      double[] x0 = vertices.get(i);
      double[] x1 = vertices.get(i == vertices.size() - 1 ? 0 : i + 1);
    
      double[] result =  new double[] { (x1[0] - x0[0]) * t + x0[0], (x1[1] - x0[1]) * t + x0[1]};
      return result;
    }
  }
}

NFlock flock;
PathMapper mapper;

FullScreen fs;
ArrayList<double[]> attractorSetup;
boolean done;

void setup() { 
  mapper = new PathMapper();
  flock = new NFlock().initialize(1280,720); 
  //fs = new FullScreen(this);
  //fs.enter();
  attractorSetup = new ArrayList<double[]>();
}

void mousePressed()
{
  if (done)
    return;
    
  mapper.addPoint(((double)mouseX) / flock.world_width, ((double)mouseY) / flock.world_height);
  fill(255);
  ellipse(mouseX - 2, mouseY - 2, 4, 4);
}

void keyReleased()
{
  if (key == 'f')
    flock.current_drive = .00000000000001;
    
  if (key == 'n')
  {
    mapper.addPath();
  }
  if (key == 'd')
  {
    mapper.done();
    
    flock.attractors = new double[mapper.paths.size()][];
    
    for (int i = 0; i < flock.attractors.length; i++)
    {
      flock.attractors[i] = new double[flock.dimension];
    }
    
    done = true;
  }
}

void keyPressed()
{
  if (key == 'f')
    flock.current_drive = flock.drive;
}

double t = 0;

void draw() { 
  if (done)
  {
    
    for (int i = 0; i < 5; i++) 
    {    
      
      for (int j = 0; j < flock.attractors.length; j++)
      {
        flock.phase += .01;
        
        if (flock.phase > Math.PI * 2)
          flock.phase = 0;
          
        t += 0.001;
        if (t > 1) t -= 1;
        
        double[] p = mapper.paths.get(j).getPosition(t);
        flock.attractors[j][0] = (p[0]);
        flock.attractors[j][1] = (p[1]);
        flock.attractors[j][2] = .5;
        flock.attractors[j][3] = .5 * Math.sin(flock.phase * flock.phase * .8 + 1) + .5;
        flock.attractors[j][4] = .5 * Math.cos(flock.phase + 1) + .5;
        flock.attractors[j][5] = .5 * Math.sin(flock.phase * .8 + 2) + .5;
      }
      flock.iterate(); 
    }
  }
}
