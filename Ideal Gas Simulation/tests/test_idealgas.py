from unittest import TestCase
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt

from idealgas import IdealGasSim, Particle

class TestIdealGas(TestCase):
    def setUp(self):
        self.sim = IdealGasSim(sim_size=100, nparticles=100, mass=1, radius=5, external_temp=100)
    
    def test_add_particles(self) -> None:
        self.sim.add_particles(10)
        
        self.assertTrue(True)
    
    def test_remove_particles(self) -> None:
        self.sim.remove_particles(10)
        
        self.assertTrue(True)
        
    def test_time_step(self) -> None:
        with self.subTest("One Time Step"):
            self.sim.time_step()
            
            self.assertTrue(True)
            
        with self.subTest("Ten Time Steps"):
            self.sim.time_step(10)   
             
            self.assertTrue(True)    
        
    def test_speed_bin(self) -> None:
        n = self.sim.speed_bin(0, 1e6)
        
        self.assertEqual(n, 100)
        
    def test_thermometer(self) -> None:
        t = self.sim.thermometer()
        
        self.assertTrue(t > 0 and t < 1000)
        
    def test_plot(self) -> None:
        fig = plt.figure(visible=False)
        
        self.sim.plot(fig)
        plt.close(fig)
        
        self.assertTrue(True)
        
    def test_check_wall_collision(self) -> None:
        with self.subTest("Particles in bounds."):
            in_bounds = True
            for particle in self.sim._particles:
                
                checks = self.sim._check_wall_collision(particle)
                if np.any(checks):
                    in_bounds = False
                    
            self.assertTrue(in_bounds)
            
        with self.subTest("Particles out of bounds."):
            for particle in self.sim._particles:
                particle.P += self.sim._sim_size
            out_of_bounds = True
            for particle in self.sim._particles:
                
                checks = self.sim._check_wall_collision(particle)
                if np.sum(checks) < 3:
                    out_of_bounds = False
                    
            self.assertTrue(out_of_bounds)
        
    def test_wall_collision(self) -> None:
        with self.subTest("No collision."):
            checks = np.asarray((False, False, False))
            particle = self.sim._particles[0]
            before = particle.V.copy()
            
            self.sim._wall_collision(particle, checks)
            after = particle.V
            
            self.assertTrue(np.all(np.equal(before, after)))
            
        with self.subTest("Collision with one wall."):
            checks = np.asarray((True, False, False))
            particle = self.sim._particles[0]
            before = particle.V[0]
            
            self.sim._wall_collision(particle, checks)
            after = particle.V[0]
            
            self.assertEqual(before, -1*after)
            
        with self.subTest("Collision with two walls."):
            checks = np.asarray((True, True, False))
            particle = self.sim._particles[0]
            before = particle.V.copy()
            
            self.sim._wall_collision(particle, checks)
            after = particle.V
            
            self.assertEqual(before[0], -1*after[0])
            self.assertEqual(before[1], -1*after[1])
            
        with self.subTest("Collision with all walls."):
            checks = np.asarray((True, True, True))
            particle = self.sim._particles[0]
            before = particle.V.copy()
            
            self.sim._wall_collision(particle, checks)
            after = particle.V
            
            self.assertTrue(np.all(np.equal(before, -1*after)))
            
    def test_temp_change(self) -> None:
        raise NotImplementedError("Test not implemented.")
        
    def test_check_particle_collision(self) -> None:
        with self.subTest("Start-up"):
            particle = self.sim._particles[0]
            result = self.sim._check_particle_collision(particle)
            self.assertIsNone(result)
            
        with self.subTest("Identical Copy"):
            particle_2 = deepcopy(particle)
            self.sim._particles.append(particle_2)
            result = self.sim._check_particle_collision(particle)
            self.assertTrue(particle_2 is result)
        
    def test_particle_collision(self) -> None:
        #Define particles 1 and 2, don't use sim
        #Run the collision!
        raise NotImplementedError("Test not implemented.")
    
    def test_reset_partition(self) -> None:
        self.sim._reset_partition()
        
        self.assertEqual(len(self.sim._partition), self.sim._n_partitions)
    
    def test_partition_index(self) -> None:
        in_partition_bounds = True
        for particle in self.sim._particles:
            
            i_x, i_y = self.sim._partition_index(particle)
            if i_x < 0 or i_x >= self.sim._n_partitions:
                in_partition_bounds = False
            if i_y < 0 or i_y >= self.sim._n_partitions:
                in_partition_bounds = False
                
        self.assertTrue(in_partition_bounds)
    
    def test_neighbor_partitions(self) -> None:
        with self.subTest("Corner"):
            N_x, N_y = self.sim._neighbor_partitions(0, 0)
            
            self.assertEqual(N_x, (0, 1))
            self.assertEqual(N_y, (0, 1))
            
        with self.subTest("Side"):
            N_x, N_y = self.sim._neighbor_partitions(0, 1)
            
            self.assertEqual(N_x, (0, 1))
            self.assertEqual(N_y, (-1, 0, 1))
            
        with self.subTest("Center"):
            N_x, N_y = self.sim._neighbor_partitions(1, 1)
            
            self.assertEqual(N_x, (-1, 0, 1))
            self.assertEqual(N_y, (-1, 0, 1))


class TestParticle(TestCase):
    def setUp(self) -> None:
        self.particle = Particle(P=np.asarray((5, 6, 7)), V=np.asarray((10, 9, 8)), m=5, r=5)
    
    def test_speed(self) -> None:
        self.assertAlmostEqual(self.particle.speed, 15.6524758)
        
    def test_KE(self) -> None:
        self.assertAlmostEqual(self.particle.KE, 612.5)
        
    def test_momentum(self) -> None:
        self.assertAlmostEqual(self.particle.momentum, 78.2623792)