function circumcenter(a, b, c) {
    const p1p2 = a.clone().sub(b);
    const p2p1 = b.clone().sub(a);
    const p1p3 = a.clone().sub(c);
    const p3p1 = c.clone().sub(a);
    const p2p3 = b.clone().sub(c);
    const p3p2 = c.clone().sub(b);

    const p1p2Xp2p3 = p1p2.clone().cross(p2p3);

    const alpha = p2p3.lengthSq() * p1p2.dot(p1p3) / (2 * p1p2Xp2p3.lengthSq());
    const beta = p1p3.lengthSq() * p2p1.dot(p2p3) / (2 * p1p2Xp2p3.lengthSq());
    const gamma = p1p2.lengthSq() * p3p1.dot(p3p2) / (2 * p1p2Xp2p3.lengthSq());

    return a.clone().multiplyScalar(alpha).add(b.clone().multiplyScalar(beta)).add(c.clone().multiplyScalar(gamma));
}

class GridCell {
    constructor(centerVertex, isAlongIcosahedronEdge) {
        this.centerVertex = centerVertex;
        this.isAlongIcosahedronEdge = isAlongIcosahedronEdge;
        this.neighbors = [];
        this.faceTriangles = [];
    }

    addNeighbor(neighbor) {
        this.neighbors.push(neighbor);
    }

    calculateFaceTriangles() {
        // NOTE: This only works correctly if the cell neigbors array has the
        // neighboring cells in a counter-clockwise order around this cell.
        const centerVertex = this.centerVertex;
        const circumcenters = this.neighbors.map(function(neighbor, i, neighbors) {
            const nextNeighbor = neighbors[(i+1) % neighbors.length];
            return circumcenter(
                centerVertex,
                neighbor.centerVertex,
                nextNeighbor.centerVertex,
            ).normalize();
        });

        this.faceTriangles = [new THREE.Triangle(
            circumcenters[0],
            circumcenters[1],
            circumcenters[2],
        )];
        for(let i = 2; i+1 < circumcenters.length; i++) {
            this.faceTriangles.push(new THREE.Triangle(
                circumcenters[0],
                circumcenters[i],
                circumcenters[i+1],
            ));
        }
    }
}

class Grid {
    constructor(N) {
        const { gridPanels, northPole, southPole } = Grid.make(N);
        this.gridPanels = gridPanels;
        this.northPole = northPole;
        this.southPole = southPole;

        const sideLength = 1 << N;
        this.size = 10*sideLength*sideLength + 2;

        // Add neighbor cells for all grid cells.
        // NOTE: it is important to add the neighbor cells in a
        // counter-clockwise direction around each cell in order for the
        // calculation of face triangles to work correctly.
        for (let i = 0; i < this.gridPanels.length; i++) {
            const gridPanel = this.gridPanels[i];
            const rightGridPanel = this.gridPanels[i < 4 ? i+1 : 0];
            const leftGridPanel = this.gridPanels[i > 0 ? i-1 : 4];

            // Handle the top corner hexagon.
            let gridCell = gridPanel[0][0];
            gridCell.addNeighbor(this.northPole);
            gridCell.addNeighbor(leftGridPanel[0][0]);
            gridCell.addNeighbor(gridPanel[1][0]);
            gridCell.addNeighbor(gridPanel[0][1]);
            gridCell.addNeighbor(rightGridPanel[1][0]);
            gridCell.addNeighbor(rightGridPanel[0][0]);

            // Add neighbor for north pole.
            // NOTE: we can do this in this loop because we loop through the
            // grid panels from left to right. We'll need to go in other
            // direction for the south pole.
            this.northPole.addNeighbor(gridCell);

            // Handle the bottom corner hexagon.
            gridCell = gridPanel[sideLength-1][2*sideLength-1];
            gridCell.addNeighbor(gridPanel[sideLength-1][2*sideLength-2]);
            gridCell.addNeighbor(leftGridPanel[sideLength-1][2*sideLength-1]);
            gridCell.addNeighbor(this.southPole);
            gridCell.addNeighbor(rightGridPanel[sideLength-1][2*sideLength-1]);
            gridCell.addNeighbor(rightGridPanel[sideLength-1][2*sideLength-2]);
            gridCell.addNeighbor(gridPanel[sideLength-2][2*sideLength-1]);

            // Handle the other 2 pentagons which are on this panel.
            gridCell = gridPanel[0][sideLength-1];
            gridCell.addNeighbor(gridPanel[0][sideLength-2]);
            gridCell.addNeighbor(gridPanel[1][sideLength-2]);
            gridCell.addNeighbor(gridPanel[1][sideLength-1]);
            gridCell.addNeighbor(gridPanel[0][sideLength]);
            gridCell.addNeighbor(rightGridPanel[sideLength-1][0]);

            gridCell = gridPanel[0][2*sideLength-1];
            gridCell.addNeighbor(gridPanel[0][2*sideLength-2]);
            gridCell.addNeighbor(gridPanel[1][2*sideLength-2]);
            gridCell.addNeighbor(gridPanel[1][2*sideLength-1]);
            gridCell.addNeighbor(rightGridPanel[sideLength-1][sideLength]);
            gridCell.addNeighbor(rightGridPanel[sideLength-1][sideLength-1]);

            // Handle the left corner.
            gridCell = gridPanel[sideLength-1][0];
            gridCell.addNeighbor(leftGridPanel[0][sideLength-2]);
            gridCell.addNeighbor(leftGridPanel[0][sideLength-1]);
            gridCell.addNeighbor(leftGridPanel[0][sideLength]);
            gridCell.addNeighbor(gridPanel[sideLength-1][1]);
            gridCell.addNeighbor(gridPanel[sideLength-2][1]);
            gridCell.addNeighbor(gridPanel[sideLength-2][0]);

            // Handle the top-left edge.
            for (let j = 1; j < sideLength-1; j++) {
                gridCell = gridPanel[j][0];
                gridCell.addNeighbor(leftGridPanel[0][j-1]);
                gridCell.addNeighbor(leftGridPanel[0][j]);
                gridCell.addNeighbor(gridPanel[j+1][0]);
                gridCell.addNeighbor(gridPanel[j][1]);
                gridCell.addNeighbor(gridPanel[j-1][1]);
                gridCell.addNeighbor(gridPanel[j-1][0]);
            }

            // Handle the top-right edge.
            for (let i = 1; i < sideLength-1; i++) {
                gridCell = gridPanel[0][i];
                gridCell.addNeighbor(gridPanel[0][i-1]);
                gridCell.addNeighbor(gridPanel[1][i-1]);
                gridCell.addNeighbor(gridPanel[1][i]);
                gridCell.addNeighbor(gridPanel[0][i+1]);
                gridCell.addNeighbor(rightGridPanel[i+1][0]);
                gridCell.addNeighbor(rightGridPanel[i][0]);
            }

            // Handle the mid-left edge.
            for (let i = 1; i < sideLength; i++) {
                gridCell = gridPanel[sideLength-1][i];
                gridCell.addNeighbor(gridPanel[sideLength-1][i-1]);
                gridCell.addNeighbor(leftGridPanel[0][sideLength+i-1]);
                gridCell.addNeighbor(leftGridPanel[0][sideLength+i]);
                gridCell.addNeighbor(gridPanel[sideLength-1][i+1]);
                gridCell.addNeighbor(gridPanel[sideLength-2][i+1]);
                gridCell.addNeighbor(gridPanel[sideLength-2][i]);
            }

            // Handle the mid-right edge.
            for (let i = 0; i < sideLength-1; i++) {
                gridCell = gridPanel[0][sideLength+i];
                gridCell.addNeighbor(gridPanel[0][sideLength+i-1]);
                gridCell.addNeighbor(gridPanel[1][sideLength+i-1]);
                gridCell.addNeighbor(gridPanel[1][sideLength+i]);
                gridCell.addNeighbor(gridPanel[0][sideLength+i+1]);
                gridCell.addNeighbor(rightGridPanel[sideLength-1][i+1]);
                gridCell.addNeighbor(rightGridPanel[sideLength-1][i]);
            }

            // Handle the bottom-left edge.
            for (let i = 0; i < sideLength-1; i++) {
                gridCell = gridPanel[sideLength-1][sideLength+i];
                gridCell.addNeighbor(gridPanel[sideLength-1][sideLength+i-1]);
                gridCell.addNeighbor(leftGridPanel[i][2*sideLength-1]);
                gridCell.addNeighbor(leftGridPanel[i+1][2*sideLength-1]);
                gridCell.addNeighbor(gridPanel[sideLength-1][sideLength+i+1]);
                gridCell.addNeighbor(gridPanel[sideLength-2][sideLength+i+1]);
                gridCell.addNeighbor(gridPanel[sideLength-2][sideLength+i]);
            }

            // Handle the bottom-right edge.
            for (let j = 1; j < sideLength-1; j++) {
                gridCell = gridPanel[j][2*sideLength-1];
                gridCell.addNeighbor(gridPanel[j][2*sideLength-2]);
                gridCell.addNeighbor(gridPanel[j+1][2*sideLength-2]);
                gridCell.addNeighbor(gridPanel[j+1][2*sideLength-1]);
                gridCell.addNeighbor(rightGridPanel[sideLength-1][sideLength+j]);
                gridCell.addNeighbor(rightGridPanel[sideLength-1][sideLength+j-1]);
                gridCell.addNeighbor(gridPanel[j-1][2*sideLength-1]);
            }

            // Handle the middle section.
            for (let j = 1; j < sideLength-1; j++) {
                for (let i = 1; i < 2*sideLength-1; i++) {
                    gridCell = gridPanel[j][i];
                    gridCell.addNeighbor(gridPanel[j][i-1]);
                    gridCell.addNeighbor(gridPanel[j+1][i-1]);
                    gridCell.addNeighbor(gridPanel[j+1][i]);
                    gridCell.addNeighbor(gridPanel[j][i+1]);
                    gridCell.addNeighbor(gridPanel[j-1][i+1]);
                    gridCell.addNeighbor(gridPanel[j-1][i]);
                }
            }
        }

        // Finally, add neighbors for the south pole.
        for (let i = 4; i >= 0; i--) {
            const gridPanel = this.gridPanels[i];
            const gridCell = gridPanel[sideLength-1][2*sideLength-1];
            this.southPole.addNeighbor(gridCell);
        }

        // Now that all of the grid cells have their neighbors all specified,
        // we can calculate the face triangles for each cell.
        for (const gridCell of this) {
            if (gridCell.neighbors.length !== 5 && gridCell.neighbors.length !== 6) {
                throw `grid cell has ${gridCell.neighbors.length} neighbors`;
            }
            gridCell.calculateFaceTriangles();
        }
    }

    calculateGridFaces() {

    }

    *[Symbol.iterator]() {
        yield this.northPole;
        const sideLength = this.gridPanels[0].length;
        for (const gridPanel of this.gridPanels) {
            for (let j = 0; j < sideLength; j++) {
                for (let i = 0; i < 2*sideLength; i++) {
                    yield gridPanel[j][i];
                }
            }
        }
        yield this.southPole;
    }

    static make(N) {
        const a = 0.525731112119133606;
        const b = 0.850650808352039932;

        // 12 Vertices of an Icosahedron.

        const vA = new THREE.Vector3(-a, 0, b);
        const vB = new THREE.Vector3( a, 0, b);
        const vC = new THREE.Vector3(-a, 0, -b);
        const vD = new THREE.Vector3( a, 0, -b);
        const vE = new THREE.Vector3( 0, b, a);
        const vF = new THREE.Vector3( 0, b, -a);
        const vG = new THREE.Vector3( 0, -b, a);
        const vH = new THREE.Vector3( 0, -b, -a);
        const vI = new THREE.Vector3( b, a, 0);
        const vJ = new THREE.Vector3(-b, a, 0);
        const vK = new THREE.Vector3( b, -a, 0);
        const vL = new THREE.Vector3(-b, -a, 0);

        const vertices = [vA, vB, vC, vD, vE, vF, vG, vH, vI, vJ, vK, vL];

        /*

        These vertices start out pretty symmetrical, each being within a plane
        of two of the x, y, or z axes. We want to rotate it about the y axis such
        that vA becomes a unit vector in the direction of the z axis. The sine of
        that rotation angle happens to be 'a' and the cosine happens to be 'b'. The
        transformation matrix for this rotation is given by:

            [ b   0  a ]
            [ 0   1  0 ]
            [ -a  0  b ]

        After applying this rotation, vA is the "north pole" and vD is the
        "south pole". This is a bit more straightforward that using a quaternion
        to rotate the vectors.

        */
        const m = new THREE.Matrix3();
        m.set( b, 0, a,
               0, 1, 0,
              -a, 0, b );

        vertices.forEach(function(vec) {
            vec.applyMatrix3(m);
        });

        // Next, we'll make the 20 triangle faces of the icosahedron.
        // The order of these vertices matters in each triangle to establish which
        // side of the triangle is the "front" side according to the right hand
        // rule. If the "back" side faces the camera then it will not be rendered
        // by default, so all of the triangles here are created such that the
        // "front" sides face out from the base icosahedron.
        // The icosahedron can be arranged into 5 logically rectangular panels,
        // each with 4 triangle faces of the icosahedron.
        const panels = [
            [vA, vB, vE, vI, vF, vD],
            [vA, vE, vJ, vF, vC, vD],
            [vA, vJ, vL, vC, vH, vD],
            [vA, vL, vG, vH, vK, vD],
            [vA, vG, vB, vK, vI, vD],
        ];

        const gridPanels = panels.map(function(panel) {
            const grid = new Array(1 << N);
            for (let i = 0; i < grid.length; i++) {
                grid[i] = new Array(1 << (N+1));
            }

            Grid.subdivideDiamond(N, 0, 0, panel, grid);

            return grid;
        });

        return {
            gridPanels: gridPanels,
            northPole: new GridCell(vA, true),
            southPole: new GridCell(vD, true),
        };
    }

    static isAlongIcosahedronEdge(sideLength, i, j) {
        return j === 0 || i === sideLength-1 || i === 2*sideLength-1 || i+j === sideLength-1 || i+j === 2*sideLength-1;
    }

    static subdivideDiamond(N, x, y, panel, grid) {
        const [a, b, c, d, e, f] = panel;

        if (N == 0) {
            grid[y][x] = new GridCell(c, Grid.isAlongIcosahedronEdge(grid.length, x, y));
            grid[y][x+1] = new GridCell(e, Grid.isAlongIcosahedronEdge(grid.length, x+1, y));
            return;
        }

        //
        //         a 
        //          /\
        //         /  \
        //      g /____\ h
        //       /\    /\
        //      /  \  /  \
        //   b /____i/____\ c
        //     \    /\    /\
        //      \  /  \  /  \
        //     j \/____k/____\ l
        //        \    /\    /\
        //         \  /  \  /  \
        //        d \/____m/____\ e
        //           \    /\    /
        //            \  /  \  /
        //           n \/____\/ o
        //              \    /
        //               \  /
        //              f \/
        //
        
        const g = a.clone().add(b).normalize();
        const h = a.clone().add(c).normalize();
        const i = b.clone().add(c).normalize();
        const j = b.clone().add(d).normalize();
        const k = c.clone().add(d).normalize();
        const l = c.clone().add(e).normalize();
        const m = d.clone().add(e).normalize();
        const n = d.clone().add(f).normalize();
        const o = e.clone().add(f).normalize();

        N--;

        let x0 = x;
        let y0 = y;
        let x1 = x + (1 << (N+1));
        let y1 = y + (1 << N);

        Grid.subdivideDiamond(N, x0, y0, [a, g, h, i, c, k], grid);
        Grid.subdivideDiamond(N, x0, y1, [g, b, i, j, k, d], grid);
        Grid.subdivideDiamond(N, x1, y0, [c, k, l, m, e, o], grid);
        Grid.subdivideDiamond(N, x1, y1, [k, d, m, n, o ,f], grid);
    }
}
