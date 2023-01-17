

@testset "Find maximum time to nil triangle" begin
    t3 = (
        CausticsEngineering.Vertex3D(5.0, 5.0, 0.0, 0.0, 0.0),
        CausticsEngineering.Vertex3D(6.0, 5.0, 0.0, 0.0, 0.0),
        CausticsEngineering.Vertex3D(5.5, 5.5, 0.0, 0.0, -3.0),
    )
    @test CausticsEngineering.find_maximum_t(t3)[1] ≈ 0.5 / 3.0 atol = 0.01


    t3 = (
        CausticsEngineering.Vertex3D(5.0, 5.0, 0.0, 0.0, 0.0),
        CausticsEngineering.Vertex3D(6.0, 5.0, 0.0, 0.0, 0.0),
        CausticsEngineering.Vertex3D(5.5, 5.5, 0.0, 0.7, -3.0),
    )
    @test CausticsEngineering.find_maximum_t(t3)[1] ≈ 0.5 / 3.0 atol = 0.01


    t3 = (
        CausticsEngineering.Vertex3D(5.0, 5.0, 0.0, 0.0, 0.0),
        CausticsEngineering.Vertex3D(6.0, 5.0, 0.0, 0.0, 0.0),
        CausticsEngineering.Vertex3D(5.5, 5.5, 0.0, -10.0, -3.0),
    )
    @test CausticsEngineering.find_maximum_t(t3)[1] ≈ 0.5 / 3.0 atol = 0.01


    t3 = (
        CausticsEngineering.Vertex3D(5.0, 5.0, 0.01, 0.1, -2.0),
        CausticsEngineering.Vertex3D(6.0, 5.0, 0.01, 0.1, 2.0),
        CausticsEngineering.Vertex3D(5.5, 5.5, 0.01, 0.5, -3.0),
    )
    @test CausticsEngineering.find_maximum_t(t3)[1] ≈ 0.1540 atol = 0.01

    t3 = (
        CausticsEngineering.Vertex3D(5.0, 5.0, 0.01, 0.1, -2.0),
        CausticsEngineering.Vertex3D(6.0, 5.0, 0.01, 0.1, 2.0),
        CausticsEngineering.Vertex3D(5.5, 5.5, 0.01, 0.5, -3.0),
    )
    @test CausticsEngineering.find_maximum_t(t3)[1] ≈ 0.1540 atol = 0.01

end
