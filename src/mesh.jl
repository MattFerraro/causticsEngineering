"""
$(SIGNATURES)

Given 3 points and their velocities, calculate the time `t` required to bring the area of that triangle to zero
"""
function find_maximum_t(p1::Vertex3D, p2::Vertex3D, p3::Vertex3D)
    # Three points A, B and C, with coordinates (x, y)
    # The area of a triangle is 1/2 * [ Ax (By - Cy) + Bx (Cy - Ay) + Cx (Ay - By)]
    # where each point of the triangle is where it will be after time t
    # i.e. a point goes from P to P+tV where V is the velocity of that point.
    # Notation is here B -> B + t_vB

    # To make the calculation simpler, everything is translated so that A is at the
    # origin of the plane and its velocity is nil.
    Br = p2.r - p1.r
    Bc = p2.c - p1.c
    Cr = p3.r - p1.r
    Cc = p3.c - p1.c

    t_vBr = p2.vr - p1.vr
    t_vBc = p2.vc - p1.vc
    t_vCr = p3.vr - p1.vr
    t_vCc = p3.vc - p1.vc

    # After this, given that Ar = Ac = t_vAr = t_vAc = 0, the area is nil iff
    # (Br + t_vBr) (Cc + t_vCc ) - (Cr + t_vCr) (Bc + t_vBc) = 0.
    # After expansion and reshuffling to have a quadratic equation where t
    # is the variable, the coefficients of that equation are:
    a = t_vCc * t_vBr - t_vBc * t_vCr
    b = -Bc * t_vCr - Cr * t_vBc + Br * t_vCc + Cc * t_vBr
    c = Br * Cc - Cr * Bc

    # if a = 0, this is just a linear equation.
    if a == 0 && b != 0
        return smallest_positive(-c / b, c / b)
    else
        discriminant = b^2 - 4a * c

        # If there is a solution
        if discriminant >= 0
            d = sqrt(discriminant)
            return smallest_positive((-b - d) / 2a, (-b + d) / 2a)
        end
    end
    # There can be no solution if, after translation, B abd C move in parallel direction.
    # C will never end up on the line AB.
    # Very unlikely with Float64.
    # Negative numbers are filtered out when calculating the minimum jiggle ratio
    return -1.0
end

find_maximum_t(p::Tuple{Vertex3D,Vertex3D,Vertex3D}) = find_maximum_t(p[1], p[2], p[3])




"""
$(SIGNATURES)

This function will take a `grid_definition x grid_definition` matrix and returns a
`grid_definition x grid_definition` mesh.

"""
function matrix_to_mesh(ϕ::Matrix{Float64})
    height, width = size(ϕ)

    mesh = FaceMesh(height, width)

    # 1 more corner than corners! Therefore compiler needs to specify exact indices.
    mesh.corners.ϕ[1:height, 1:width] .= ϕ[1:height, 1:width]

    # The borders' height is forced at 0.
    reset_border_values!(mesh.corners)

    return mesh
end


"""
$(SIGNATURES)

Create a mash of the solid block to be carved.
Imagining a mid-plane where the heightmap is carved, start with a block whose lower plane is
bottom_distance away and the highr plane is top_distance above.
The distances are measured in meters and converted to pixels when necessary.

The final mesh is a vector of vertices each with an ID number and a position.

"""
function create_solid_as_dict(
    mesh::FaceMesh;
    bottom_distance = Bottom_Offset,
    top_distance = Top_Offset,
)

    height, width = size(mesh)
    bottom_offset = bottom_distance / Meters_Per_Pixel
    top_offset = top_distance / Meters_Per_Pixel

    # To do this, we use several dictionaries before saving
    # Dictionary 1: Triangle to 3 strings each naming a triangle vertex.
    # Dictionary 2: Named vertex to its x, y, z coordinates.
    # Dictionary 3: Named vertex to a unique numerical index.
    triangle_dict = Dict{String,Vector{String}}()
    vertex_dict = Dict{String,Tuple{Float64,Float64,Float64}}()

    mesh_height, mesh_width = size(mesh.corners.ϕ)
    for row ∈ 1:mesh_height-1, col ∈ 1:mesh_width-1

        # Top face, add vertices, top and bottom triangles
        for side ∈ [:top, :bottom]
            side_name = side == :top ? "Top" : "Bottom"
            t = side == :top ? mesh.toptriangles[row, col] : mesh.bottriangles[row, col]

            # Add the vertices 1 by 1
            vertices_names = ["", "", ""]
            for i ∈ 1:3
                v_row, v_col = t[i]
                r = mesh.corners.r[v_row, v_col]
                c = mesh.corners.c[v_row, v_col]
                ϕ = mesh.corners.ϕ[v_row, v_col]
                name = vertex_name("Top", row, col)
                vertices_names[i] = name
                vertex_dict[name] = (r, c, ϕ)
            end

            # Add the triangle
            triangle_dict[triangle_name(side_name, row, col, side)] = vertices_names
        end

        # Bottom face
        v_coord_top = [(row, col), (row + 1, col), (row, col + 1)]
        vertices_names_top = ["", "", ""]

        v_coord_bot = [(row, col + 1), (row + 1, col), (row + 1, col + 1)]
        vertices_names_bot = ["", "", ""]

        # Add the vertices 1 by one
        for i ∈ 1:3
            (v_row, v_col) = v_coord_top[i]
            name = vertex_name("Bottom", v_row, v_col)
            vertices_names_top[i] = name
            vertex_dict[name] =
                (Float64(v_row), Float64(v_col), Float64(-bottom_offset * 0))

            (v_row, v_col) = v_coord_bot[i]
            name = vertex_name("Bottom", v_row, v_col)
            vertices_names_bot[i] = name
            vertex_dict[name] = (Float64(v_row), Float64(v_col), Float64(-bottom_offset))
        end

        # Add the triangles
        triangle_dict[triangle_name("Bottom", row, col, :top)] = vertices_names_top
        triangle_dict[triangle_name("Bottom", row, col, :bottom)] = vertices_names_bot
    end


    ###
    ### Side meshes
    ###
    # Build triangles to create side meshes and close the mesh
    # The mesh is made of a single couple top/bottom triangles joining the bottom face
    # and the top face.
    # All the vertices already exist

    # (Looking from above), west/east sides.
    for (col, face_name) ∈ [(1, "West"), (width, "East")]
        for row = 1:height-1

            v_coord_top = [(:top, row, col), (:bottom, row, col), (:top, row + 1, col)]
            vertices_names_top = ["", "", ""]
            for i ∈ 1:3
                vertices_names_top[i] = vertex_name(
                    (v_coord_top[i][1] == :top ? "Top" : "Bottom"),
                    v_coord_top[i][2],
                    v_coord_top[i][3],
                )
            end

            v_coord_bot =
                [(:top, row + 1, col), (:bottom, row + 1, col), (:bottom, row, col)]
            vertices_names_bot = ["", "", ""]
            for i ∈ 1:3
                vertices_names_bot[i] = vertex_name(
                    (v_coord_bot[i][1] == :top ? "Top" : "Bottom"),
                    v_coord_bot[i][2],
                    v_coord_bot[i][3],
                )
            end

            # Add the triangles
            triangle_dict[triangle_name(face_name, row, col, :top)] = vertices_names_top
            triangle_dict[triangle_name(face_name, row, col, :bottom)] = vertices_names_bot
        end
    end

    # (Looking from above), north/south sides.
    for (row, face_name) ∈ [(1, "North"), (height, "South")]
        for col ∈ 1:width-1

            v_coord_top = [(:top, row, col), (:bottom, row, col), (:top, row, col + 1)]
            vertices_names_top = ["", "", ""]
            for i ∈ 1:3
                vertices_names_top[i] = vertex_name(
                    (v_coord_top[i][1] == :top ? "Top" : "Bottom"),
                    v_coord_top[i][2],
                    v_coord_top[i][3],
                )
            end

            v_coord_bot = [(:top, row, col + 1), (:bottom, row, col + 1), (:top, row, col)]
            vertices_names_bot = ["", "", ""]
            for i ∈ 1:3
                vertices_names_bot[i] = vertex_name(
                    (v_coord_bot[i][1] == :top ? "Top" : "Bottom"),
                    v_coord_bot[i][2],
                    v_coord_bot[i][3],
                )
            end


            # Add the triangles
            triangle_dict[triangle_name(face_name, row, col, :top)] = vertices_names_top
            triangle_dict[triangle_name(face_name, row, col, :bottom)] = vertices_names_bot
        end
    end

    return triangle_dict, vertex_dict
end


"""
$(SIGNATURES)

"""
function create_solid(
    mesh::FaceMesh;
    bottom_distance = Bottom_Offset,
    top_distance = Top_Offset,
)

    triangle_dict, vertex_dict = create_solid_as_dict(
        mesh;
        bottom_distance = Bottom_Offset,
        top_distance = Top_Offset,
    )

    _, _, triangle_index, vertex_index = create_solid(
        triangle_dict,
        vertex_dict;
        bottom_distance = bottom_distance,
        top_distance = top_distance,
    )

    return triangle_dict, vertex_dict, triangle_index, vertex_index
end



"""
$(SIGNATURES)

"""
function create_solid(
    triangle_dict::Dict{String,Vector{String}},
    vertex_dict::Dict{String,Tuple{Float64,Float64,Float64}};
    bottom_distance = Bottom_Offset,
    top_distance = Top_Offset,
)

    index_dict = Dict{String,Int64}()
    index_of_vertex = 0
    for name ∈ keys(vertex_dict)
        index_of_vertex += 1
        index_dict[name] = index_of_vertex
    end

    vertex_index = Vector{Tuple{Float64,Float64,Float64}}(undef, length(vertex_dict))
    # For each vertex using its name and coordinates
    for (name, v) ∈ vertex_dict
        # Determine the index of the vertex
        i = index_dict[name]
        x = v[1]
        y = v[2]
        z = v[3]
        vertex_index[i] = (x, y, z)
    end

    # The connections are listed as a tuple of 3 arrays.
    # First array contains the first point of a triangle, and so on
    triangle_index = [
        (index_dict[v[1]], index_dict[v[2]], index_dict[v[3]]) for
        v ∈ values(triangle_dict)
    ]

    return triangle_dict, vertex_dict, triangle_index, vertex_index

end
