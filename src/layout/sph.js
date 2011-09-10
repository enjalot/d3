// A rudimentary sph layout
d3.layout.sph = function() {
  var force = {},
      event = d3.dispatch("tick"),
      size = [1, 1],
      drag,
      alpha,
      friction = .0,
      linkDistance = d3_layout_sphLinkDistance,
      linkStrength = d3_layout_sphLinkStrength,

      //constants
      //from force.js
      charge = -30,
      gravity = .1,
      theta = .8,
      //sph
      rho0 = 1000,                  //rest density [ kg/m^3 ]
      VF = .0262144,                //simulation volume [ m^3 ]
      max_num = 4096,               //won't actually simulate with this many
      interval,
      nodes = [],
      links = [],
      distances,
      test,
      strengths;

    //initialize from constants
    var VP = VF / max_num;              //particle volume [ m^3 ]
    var m = rho0 * VP;                  //particle mass [ kg ]
    var re = Math.pow(VP,(1/3));               //particle radius [ m ]
    var rest_distance = .87 * re;    //rest distance between particles [ m ]


    //hardcoding the simulation domain for now
    var min = { x: 0, y: 0 }; 
    var max = { x: 50, y: 50 }; 
    var V = 50 * 50 * 50;       //Volume

    var smoothing_radius = 2.0 * rest_distance;      //smoothing radius for SPH Kernels
    var boundary_distance = .5 * rest_distance;      //for calculating collision with boundary
    var sim_scale = Math.pow((VF / V),(1/3));         //[m^3 / world m^3 ]

    /*
    sx = max.x - min.x; 
    sy = max.y - min.y; 
    cell_size = smoothing_radius / sim_scale;
    resx = Math.ceil(sx / cell_size);
    resy = Math.ceil(sy / cell_size);
    sx = resx * cell_size;
    sy = resy * cell_size;
    max.x = min.x + sx; 
    max.y = min.y + sy; 
    //nb_cells = int(resx *resy);
    //console.log(nb_cells);
    */
 

    // Find the nodes within cell
    function find(quadtree, x0, y0, x3, y3) {
      var points = [];
      quadtree.visit(function(node, x1, y1, x2, y2) {
        var p = node.point;
        if (p && (p.x >= x0) && (p.x < x3) && (p.y >= y0) && (p.y < y3)) points.push(p);
        return x1 >= x3 || y1 >= y3 || x2 < x0 || y2 < y0;
      });
      return points;
    }



  function repulse(node, kc) {
    return function(quad, x1, y1, x2, y2) {
      if (quad.point !== node) {
        var dx = quad.cx - node.x,
            dy = quad.cy - node.y,
            dn = 1 / Math.sqrt(dx * dx + dy * dy);
            if( node == 0)
            {
                nodes[node].color = 1;
                console.log("quad " + quad + " " + 1 / dn);
            }

        /* Barnes-Hut criterion. */
        if ((x2 - x1) * dn < theta) {
          var k = kc * quad.count * dn * dn;
          node.px -= dx * k;
          node.py -= dy * k;
          return true;
        }

        if (quad.point && isFinite(dn)) {
          var k = kc * dn * dn;
          node.px -= dx * k;
          node.py -= dy * k;
        }
      }
    };
  }

  function tick() {
    var n = nodes.length,
        m = links.length,
        q,
        i, // current index
        o, // current object
        s, // current source
        t, // current target
        l, // current distance
        k, // current force
        x, // x-distance
        y; // y-distance

    // gauss-seidel relaxation for links
    // brute force density
/*
    for (i = 0; i < m; ++i) {
      o = links[i];
      s = o.source;
      t = o.target;
      x = t.x - s.x;
      y = t.y - s.y;
      if (l = (x * x + y * y)) {
        l = alpha * strengths[i] * ((l = Math.sqrt(l)) - distances[i]) / l;
        x *= l;
        y *= l;
        t.x -= x * (k = s.weight / (t.weight + s.weight));
        t.y -= y * k;
        s.x += x * (k = 1 - k);
        s.y += y * k;
      }
    }
*/
    // apply gravity forces
/*
    if (k = alpha * gravity) {
      x = size[0] / 2;
      y = size[1] / 2;
      //y = size[1];
      i = -1; if (k) while (++i < n) {
        o = nodes[i];
        o.x += (x - o.x) * k;
        o.y += (y - o.y) * k;
        //o.y -= k
      }
    }
*/

    // compute quadtree center of mass and apply charge forces
/*
    if (k = alpha * charge) {
      d3_layout_forceAccumulate(q = d3.geom.quadtree(nodes));
      i = -1; while (++i < n) {
        if (!(o = nodes[i]).fixed) {
          o.color = 0;
          q.visit(repulse(o, k));
        }
      }
    }
*/
    nodes.forEach(function(n) { n.color = 0; } );
    n0 = nodes[0];
    n0.color = 1;
    console.log("n0");
    console.log(n0);
    qt = d3.geom.quadtree(nodes);
    nn = find(qt, n0.x - 50, n0.y - 50, n0.x + 50, n0.y + 50);
    console.log("nn");
    console.log(nn);
    nn.forEach(function(d) { if (d != n0) { d.color = 2; }});

    // position verlet integration
//
    i = -1; while (++i < n) {
      o = nodes[i];
      //console.log("verlet: " + o.color);
      if (o.fixed) {
        o.x = o.px;
        o.y = o.py;
      } else {
        o.x -= (o.px - (o.px = o.x)) * friction;
        o.y -= (o.py - (o.py = o.y)) * friction;
      }
    }
//

    event.tick.dispatch({type: "tick", alpha: alpha});

    //alpha *= .1;
    console.log(alpha);
    //console.log( (alpha *= .99) < .005);
    // simulated annealing, basically
    //return (alpha *= .99) < .005;
    //return alpha < .001
    return true;
    //return false;
  }

  force.on = function(type, listener) {
    event[type].add(listener);
    return force;
  };

  force.nodes = function(x) {
    if (!arguments.length) return nodes;
    nodes = x;
    return force;
  };

  force.links = function(x) {
    if (!arguments.length) return links;
    links = x;
    return force;
  };

  force.size = function(x) {
    if (!arguments.length) return size;
    size = x;
    return force;
  };

  force.linkDistance = function(x) {
    if (!arguments.length) return linkDistance;
    linkDistance = d3.functor(x);
    return force;
  };

  // For backwards-compatibility.
  force.distance = force.linkDistance;

  force.linkStrength = function(x) {
    if (!arguments.length) return linkStrength;
    linkStrength = d3.functor(x);
    return force;
  };

  force.friction = function(x) {
    if (!arguments.length) return friction;
    friction = x;
    return force;
  };

  force.charge = function(x) {
    if (!arguments.length) return charge;
    charge = x;
    return force;
  };

  force.gravity = function(x) {
    if (!arguments.length) return gravity;
    gravity = x;
    return force;
  };

  force.theta = function(x) {
    if (!arguments.length) return theta;
    theta = x;
    return force;
  };

  force.start = function() {
    var i,
        j,
        n = nodes.length,
        m = links.length,
        w = size[0],
        h = size[1],
        neighbors,
        o;

    for (i = 0; i < n; ++i) {
      (o = nodes[i]).index = i;
      o.weight = 0;
      o.color = 0;
    }
test = 10;

    distances = [];
    strengths = [];
    for (i = 0; i < m; ++i) {
      o = links[i];
      if (typeof o.source == "number") o.source = nodes[o.source];
      if (typeof o.target == "number") o.target = nodes[o.target];
      distances[i] = linkDistance.call(this, o, i);
      strengths[i] = linkStrength.call(this, o, i);
      ++o.source.weight;
      ++o.target.weight;
    }

    for (i = 0; i < n; ++i) {
      o = nodes[i];
      if (isNaN(o.x)) o.x = position("x", w);
      if (isNaN(o.y)) o.y = position("y", h);
      if (isNaN(o.px)) o.px = o.x;
      if (isNaN(o.py)) o.py = o.y;
    }

    // initialize node position based on first neighbor
    function position(dimension, size) {
      var neighbors = neighbor(i),
          j = -1,
          m = neighbors.length,
          x;
      while (++j < m) if (!isNaN(x = neighbors[j][dimension])) return x;
      return Math.random() * size;
    }

    // initialize neighbors lazily
    function neighbor() {
      if (!neighbors) {
        neighbors = [];
        for (j = 0; j < n; ++j) {
          neighbors[j] = [];
        }
        for (j = 0; j < m; ++j) {
          var o = links[j];
          neighbors[o.source.index].push(o.target);
          neighbors[o.target.index].push(o.source);
        }
      }
      return neighbors[i];
    }

    return force.resume();
  };

  force.resume = function() {
    alpha = .001;
    d3.timer(tick);
    return force;
  };

  force.stop = function() {
    alpha = 0;
    return force;
  };

  force.getAlpha = function() {
    return alpha;
  };
  

  



  // use `node.call(force.drag)` to make nodes draggable
  force.drag = function() {
    if (!drag) drag = d3.behavior.drag()
        .on("dragstart", dragstart)
        .on("drag", d3_layout_forceDrag)
        .on("dragend", d3_layout_forceDragEnd);

    this.on("mouseover.force", d3_layout_forceDragOver)
        .on("mouseout.force", d3_layout_forceDragOut)
        .call(drag);
  };

  function dragstart(d) {
    d3_layout_forceDragOver(d3_layout_forceDragNode = d);
    d3_layout_forceDragForce = force;
  }

  return force;
};

var d3_layout_forceDragForce,
    d3_layout_forceDragNode;

function d3_layout_forceDragOver(d) {
  d.fixed |= 2;
}

function d3_layout_forceDragOut(d) {
  if (d !== d3_layout_forceDragNode) d.fixed &= 1;
}

function d3_layout_forceDragEnd() {
  d3_layout_forceDrag();
  d3_layout_forceDragNode.fixed &= 1;
  d3_layout_forceDragForce = d3_layout_forceDragNode = null;
}

function d3_layout_forceDrag() {
  d3_layout_forceDragNode.px += d3.event.dx;
  d3_layout_forceDragNode.py += d3.event.dy;
  d3_layout_forceDragForce.resume(); // restart annealing
}

function d3_layout_forceAccumulate(quad) {
  var cx = 0,
      cy = 0;
  quad.count = 0;
  if (!quad.leaf) {
    var nodes = quad.nodes,
        n = nodes.length,
        i = -1,
        c;
    while (++i < n) {
      c = nodes[i];
      if (c == null) continue;
      d3_layout_forceAccumulate(c);
      quad.count += c.count;
      cx += c.count * c.cx;
      cy += c.count * c.cy;
    }
  }
  if (quad.point) {
    // jitter internal nodes that are coincident
    if (!quad.leaf) {
      quad.point.x += Math.random() - .5;
      quad.point.y += Math.random() - .5;
    }
    quad.count++;
    cx += quad.point.x;
    cy += quad.point.y;
  }
  quad.cx = cx / quad.count;
  quad.cy = cy / quad.count;
}

function d3_layout_sphLinkDistance(link) {
  return 20;
}

function d3_layout_sphLinkStrength(link) {
  return 1;
}
