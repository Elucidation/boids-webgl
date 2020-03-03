

/*================Creating a canvas=================*/
var canvas = document.getElementById('my_Canvas');
var toggleButton = document.getElementById('toggle');
// gl = canvas.getContext('experimental-webgl',  {antialias: false}); 
gl = canvas.getContext('experimental-webgl'); 

/*==========Defining and storing the geometry=======*/
var N_BODS = 1000;

document.getElementById('num_objs').innerText = N_BODS;

var renderTimeElement = document.getElementById('render_time');
var stepTimeElement = document.getElementById('step_time');
var comparisonsElement = document.getElementById('comparisons');
var bruteComparisonsElement = document.getElementById('brute_comparisons');
var reductionElement = document.getElementById('reduction');

bruteComparisonsElement.innerText = N_BODS*N_BODS;

// The preferred dimensions inside canvas.
var bounds = {x:-2.0, y:-1.0, x1:2.0, y1:1.0};
var grid_size = {x:50,y:50}; // divide space up into grid cells.

// interleaved x y
let positions = new Float32Array(N_BODS*2);
// let velocities = new Float32Array(N_BODS*2);
// interleaved dx, dy, d2
let distances = new Float32Array((N_BODS*N_BODS)*3);

// Bodies struct containing all bodies
var bods;

function resetBodies() {
  if (bods) {
    bods.pos = null;
    bods.vel = null;
    bods.mass = null;
  }

  bods = {
    pos:{x:new Float32Array(N_BODS),y:new Float32Array(N_BODS)},
    vel:{x:new Float32Array(N_BODS),y:new Float32Array(N_BODS)},
    N:N_BODS
  };
}

function randFloat(lower, upper) {
  return Math.random() * (upper - lower) + lower;
}


function initPoints() {
  let vmax = 0.5;
  for (var i = 0; i < N_BODS; i++) {
    bods.pos.x [i] = randFloat(bounds.x, bounds.x1) ; // x;
    bods.pos.y [i] = randFloat(bounds.y, bounds.y1); // y;
    bods.vel.x [i] = randFloat(-vmax, vmax); // vx
    bods.vel.y [i] = randFloat(-vmax, vmax); // vy

    positions[i*2]   = bods.pos.x[i];
    positions[i*2+1] = bods.pos.y[i];
  }
}

resetBodies();
initPoints();

// Create an empty buffer object to store the vertex buffer
var vertex_buffer = gl.createBuffer();

// Set the view port
gl.viewport(0,0,canvas.width,canvas.height);


/*=========================Shaders========================*/
loadShaders();

function loadShaders() {
  let vertCode = `
  attribute vec2 coordinates;
  uniform mat4 u_xformMatrix;

  void main(void) {
    gl_Position = u_xformMatrix * vec4(coordinates, 0.0, 0.5);
    gl_PointSize = 1.0;
  }
  `;

  let vertShader = gl.createShader(gl.VERTEX_SHADER);
  gl.shaderSource(vertShader, vertCode);
  gl.compileShader(vertShader);

  // fragment shader source code
  let fragCode = `
  void main(void) {
    gl_FragColor = vec4(1.0, 1.0, 1.0, 0.5);
  }
  `;

  let fragShader = gl.createShader(gl.FRAGMENT_SHADER);
  gl.shaderSource(fragShader, fragCode);
  gl.compileShader(fragShader);

  let shaderProgram = gl.createProgram();
  gl.attachShader(shaderProgram, vertShader); 
  gl.attachShader(shaderProgram, fragShader);
  gl.linkProgram(shaderProgram);
  gl.useProgram(shaderProgram);

  /*===================scaling==========================*/

  var Sx =0.25; 
  var Sy = Sz = 0.5;
  var xformMatrix = new Float32Array([
    Sx,   0.0,  0.0,  0.0,
    0.0,  Sy,   0.0,  0.0,
    0.0,  0.0,  Sz,   0.0,
    0.0,  0.0,  0.0,  1.0  
  ]);

  var u_xformMatrix = gl.getUniformLocation(shaderProgram, 'u_xformMatrix');
  gl.uniformMatrix4fv(u_xformMatrix, false, xformMatrix);

  /*======== Associating shaders to buffer objects ========*/

  // Bind vertex buffer object
  gl.bindBuffer(gl.ARRAY_BUFFER, vertex_buffer);

  // Set attribute coordinates to have stride 2 etc.
  var coord = gl.getAttribLocation(shaderProgram, "coordinates");
  // gl.vertexAttribPointer(attrib, index, gl.FLOAT, false, stride, offset);
  gl.vertexAttribPointer(coord, 2, gl.FLOAT, false, 0, 0);
  // Enable the attribute
  gl.enableVertexAttribArray(coord);
}



/*============= Drawing ===============*/
var renderTime;
function render() {
  let startTime = window.performance.now();
  // Push updated positions to GPU.
  // Bind appropriate array buffer to it
  gl.bindBuffer(gl.ARRAY_BUFFER, vertex_buffer);
  // Pass the vertex data to the buffer
  gl.bufferData(gl.ARRAY_BUFFER, positions, gl.DYNAMIC_DRAW);
  // Unbind the buffer
  gl.bindBuffer(gl.ARRAY_BUFFER, null);

  // Clear the canvas
  gl.clearColor(0.0, 0.0, 0.0, 1.0);
  // Clear the color buffer bit
  gl.clear(gl.COLOR_BUFFER_BIT);

  // Draw the points
  gl.drawArrays(gl.POINTS, 0, N_BODS);

  renderTime  = window.performance.now() - startTime;
  renderTimeElement.innerText = (renderTime).toFixed(3);
}

render(); // Render once to show something on screen.


/*==== Update state logic ===*/

// Around 5-10ms for 1000 objects.
function updateDistances() {
  // updates distances as dx, dy, d^2 as a interleaved float array of length N*3
  for (var i = bods.N - 1; i >= 1; i--) {
    for (var j = i - 1; j >= 0; j--) {
      let dx = bods.pos.x[j] - bods.pos.x[i];
      let dy = bods.pos.y[j] - bods.pos.y[i];
      let d2 = dx*dx + dy*dy;
      // let dx = positions[2*j] - positions[2*i];
      // let dy = positions[2*j+1] - positions[2*i+1];
      // let d2 = dx*dx + dy*dy;

      let didx = (bods.N*i + j)*3;
      distances[didx] = dx;
      distances[didx+1] = dy;
      distances[didx+2] = d2;

      // reverse direction distance.
      let didx2 = (bods.N*j + i)*3;
      distances[didx2] = -dx;
      distances[didx2+1] = -dy;
      distances[didx2+2] = d2;
    }
  }
}

function getCell(x, y, bounds, gridSize) {
  let cx = Math.floor(((x - bounds.x) / (bounds.x1 - bounds.x)) * gridSize.x);
  let cy = Math.floor(((y - bounds.y) / (bounds.y1 - bounds.y)) * gridSize.y);
  cx = cx >= gridSize.x ? gridSize.x-1 : cx;
  cy = cy >= gridSize.y ? gridSize.y-1 : cy;
  return {x: cx, y: cy};  
}

// grid = buildSpatialGrid(bods.pos, bounds, grid_size);
function buildSpatialGrid(positions, bounds, gridSize) {
  // Given a set of positions {x,y}, bounds {x,y,x1,y1} and gridSize {x,y}, 
  // return a 2D array of size [x,y] with each cell containing a list of position indices
  // corresponding to the position index that is in that cell.
  // example: grid = buildSpatialGrid(bods.pos, 
  //                                  {x:-1.0, y:-1.0, x1:1.0, y1:1.0},
  //                                  {x:10,y:10})
  // Init grid with gridSize.y rows and gridSize.x columns.
  let grid = Array(gridSize.y).fill(0).map(x => Array(gridSize.x).fill(0).map(x => Array()));

  // Now, for each position, insert into appropriate grid cell.
  for (var i = 0; i < positions.x.length; i++) {
    let cell = getCell(positions.x[i], positions.y[i], bounds, gridSize);
    grid[cell.y][cell.x].push(i);
  }

  return grid;
}

function getDistSquared(i, j, positions) {
  let dx = positions[2*i] - positions[2*j];
  let dy = positions[2*i+1] - positions[2*j+1];
  return dx*dx + dy*dy
}

function getCentroid(bodsvec, i, neighbors) {
  let sx = bodsvec.x[i];
  let sy = bodsvec.y[i];
  for (var k = 0; k < neighbors.length; k++) {
    let oi = neighbors[k];
    sx += bodsvec.x[oi];
    sy += bodsvec.y[oi];
  }
  // return average position of all neighbors and i.
  return {x:sx/(neighbors.length+1), y:sy/(neighbors.length+1)};
}

var grid;

function boidsUpdate() {
  const neighborRadius = 0.04**2;
  const personalSpaceRadius = 0.02**2;

  const oldWeight = 1.0;
  const centroidWeight = 0.01;
  const personalSpaceWeight = 0.0001;
  const headingMatchWeight = 0.1;
  const randWeight = 0.2;

  const neighborOffsets = [[0,0],[-1,-1],[-1,0],[-1,1],[0,1],[1,1],[1,0],[1,-1],[0,-1]];
  
  let compareCount = 0;

  // updateDistances();
  grid = buildSpatialGrid(bods.pos, bounds, grid_size);
  
  let posx = bods.pos.x;
  let posy = bods.pos.y;
  let velx = bods.vel.x;
  let vely = bods.vel.y;

  for (var i = 0; i < bods.N; i++) {
    // Boids
    let neighbors = [];
    // let tooCloseNeighbors = [];
    let cellNeighbors = [];

    // Find all cell neighbors using spatial grid (includes self).
    let ctrCellIdx = getCell(bods.pos.x[i], bods.pos.y[i], bounds, grid_size);
    for (var ni = 0; ni < neighborOffsets.length; ni++) {
      let neighborCell = {x:ctrCellIdx.x + neighborOffsets[ni][0], 
                          y:ctrCellIdx.y + neighborOffsets[ni][1]};
      if (neighborCell.x >=0 && neighborCell.x < grid_size.x &&
          neighborCell.y >=0 && neighborCell.y < grid_size.y) {
        // valid neighbor cell add to neighbors to check.
        let cell = grid[neighborCell.y][neighborCell.x];
        cellNeighbors = cellNeighbors.concat(cell);
      }
    }

    // Move away from boids too close.
    let closeVec = {x: 0, y: 0};

    // Iterate only over bodies in cell neighbors.
    

    for (var cidx = 0; cidx < cellNeighbors.length; cidx++) {
      j = cellNeighbors[cidx];
      if (i == j) {
        continue;
      }
      compareCount += 1;
      // let d2 = getDistSquared(j, i, positions);
      // let didx = (bods.N*i + j)*3;
      // let d2 = distances[didx+2];

      let dx = posx[j] - posx[i];
      let dy = posy[j] - posy[i];
      let d2 = dx*dx + dy*dy;


      if (d2 < neighborRadius) {
        neighbors.push(j);
        if (d2 < personalSpaceRadius) {
          // tooCloseNeighbors.push(j);
          closeVec.x -= dx / d2;
          closeVec.y -= dy / d2;
        }
      } 
    }

    // Move towards centroid of boids in neighborhood.
    let centroid = getCentroid(bods.pos, i, neighbors);
    let centroidVec = {x: (centroid.x - posx[i]),
                       y: (centroid.y - posy[i])};

    // Match velocity heading to boids in neighborhood.
    let headingVec = getCentroid(bods.vel, i, neighbors);    

    
    let randVec = {x: randFloat(-1.0, 1.0), y: randFloat(-1.0, 1.0)};

    let vx = oldWeight*bods.vel.x[i] + centroidVec.x*centroidWeight + closeVec.x*personalSpaceWeight + headingVec.x * headingMatchWeight + randVec.x*randWeight;
    let vy = oldWeight*bods.vel.y[i] + centroidVec.y*centroidWeight + closeVec.y*personalSpaceWeight + headingVec.y * headingMatchWeight + randVec.y*randWeight;
    let vnorm = Math.sqrt(vx*vx+vy*vy);
    velx[i] = vx / vnorm;
    vely[i] = vy / vnorm;
  }

  // console.log(compareCount);
  comparisonsElement.innerText = compareCount;
  reductionElement.innerText = (100.0* compareCount / (N_BODS*N_BODS)).toFixed(2);
}

function updatePositions(dt) {
  // for (var i = 0; i < N_BODS; i++) {
  //   positions[i*2] += velocities[i*2] * dt;
  //   positions[i*2+1] += velocities[i*2+1] * dt;
  // }
  for (var i=0;i<bods.N;i++) {
    bods.pos.x[i] += bods.vel.x[i]*dt;
    bods.pos.y[i] += bods.vel.y[i]*dt;

    positions[i*2] = bods.pos.x[i];
    positions[i*2+1] = bods.pos.y[i];
  }
}

function checkBoundaryConditions(bods, bounds) {
  let bouncePower = .9; // friction on bouncing off of bounds.
  for (var i = 0; i < bods.N; i++) {
    // Bounce off walls.
    // X
    if (bods.pos.x[i] < bounds.x) {
      bods.pos.x[i] = bounds.x;
      bods.vel.x[i] *= -bouncePower;
    } 
    else if (bods.pos.x[i] > bounds.x1) {
      bods.pos.x[i] = bounds.x1;
      bods.vel.x[i] *= -bouncePower;
    }
    // Y
    if (bods.pos.y[i] < bounds.y) {
      bods.pos.y[i] = bounds.y;
      bods.vel.y[i] *= -bouncePower;
    } 
    else if (bods.pos.y[i] > bounds.y1) {
      bods.pos.y[i] = bounds.y1;
      bods.vel.y[i] *= -bouncePower;
    }
  }
}

var stepTime = 0;
function step() {
  let startTime = window.performance.now();
  boidsUpdate();
  updatePositions(timeDelta);

  checkBoundaryConditions(bods, bounds);

  stepTime  = stepTime*0.9 + 0.1*(window.performance.now() - startTime);
  stepTimeElement.innerText = stepTime.toFixed(3);
}

var isRunning = false;
var drawInterval = null;
var physicsInterval = null;
var timeDelta = 0.003;

toggleRun(); // Toggle to start running by default.

function toggleRun() {
  isRunning = !isRunning;
  if (isRunning) {
    drawInterval = setInterval(render, 20); // 25 fps
    physicsInterval = setInterval(step, 20); // 10 fps
    console.log("Start");
    toggleButton.innerText = "Running";
  } else {
    clearInterval(drawInterval);
    clearInterval(physicsInterval);
    console.log("Stop");
    toggleButton.innerText = "Stopped";
  }
}
