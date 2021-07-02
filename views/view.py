# CENG 488 Assignment8 by
# Elman Hamdi
# 240201036
# June 2021
import copy

from utils import *
from objects import *
import numpy as np
from multiprocessing import Pool
from functools import partial

NOF_REFRACTION = 5
NOF_REFLECTION = 5
NOF_MAX_BOUNCE = 5

class View:
    def __init__(self, camera, scene=None, bgColor=Color(0, 0, 0, 1), subsample = 1, ambient_sample = 5, max_bounce = 5, num_refle = 5, num_refra = 5):
        self.camera = camera
        self.scene = scene
        self.bgColor = bgColor
        self.subsample = 1 if subsample < 1 else int(subsample)
        self.ambient_sample = 1 if ambient_sample < 1 else int(ambient_sample)

        #global NOF_MAX_BOUNCE
        #NOF_MAX_BOUNCE = 1 if max_bounce < 1 else int(max_bounce)

        global NOF_REFRACTION
        NOF_REFRACTION = 1 if num_refra < 1 else int(num_refra)

        global NOF_REFLECTION
        NOF_REFLECTION = 1 if num_refle < 1 else int(num_refle)


    def calculate_ray_tracing_with_thread_on_qt(self, starting_line, nof_thread=5):

        func = partial(self.rt_for_qt_thread, nof_thread)

        lst_tmp= []
        for i in range(nof_thread):
            lst_tmp.append(starting_line + i)

        with Pool(processes=nof_thread) as pool:
            results = pool.map(func, lst_tmp)

        lst = results[0]
        for i in range(1, len(results)):
            lst = lst + results[i]
        return lst

    def rt_for_qt_thread(self,  num_thread=5,starting_line = 1):
        cam_h = self.camera.window.height
        cam_w = self.camera.window.width
        pixel_list = np.full((num_thread, cam_w, 3), fill_value=self.bgColor.getRGB(), dtype=np.uint8)

        cam_front = Vec3d.positions_to_vec(self.camera.center, self.camera.eye, ).normalize()
        window_right_vec = Vec3d.cross_product(self.camera.up, cam_front).normalize()
        window_up_vec = Vec3d.cross_product(cam_front, window_right_vec, ).normalize()

        y = starting_line
        for x in range(cam_w):
            r = window_right_vec
            u = window_up_vec
            pos_on_window = ((r * (-cam_w / 2 + x)).add(u * (cam_h / 2 - y))).add(
                cam_front * self.camera.window_distance)
            pixel_pos_in_space = Pos3d.add(self.camera.eye, pos_on_window)

            ray_direction = HomogeneusCoor.normalize(
                Vec3d.positions_to_vec(pixel_pos_in_space, self.camera.eye))
            ray = Ray(start_pos=pixel_pos_in_space, direction=ray_direction)

            # Apply Antialising
            color = Color(r=0, g=0, b=0)
            if self.subsample != 0:
                starting_point = pos_on_window.sub(window_right_vec * (1 / 2)).add(window_up_vec * (1 / 2))
                step = int(self.subsample ** (1 / 2))
                dis_r = window_right_vec * (1 / step)
                dis_u = window_up_vec * (-1 / step)
                for s1 in range(1, step + 1):
                    for s2 in range(1, step + 1):
                        new_pos_on_window = starting_point.add(dis_r * s1).add(dis_u * s2)
                        pixel_pos_in_space = Pos3d.add(self.camera.eye, new_pos_on_window)
                        ray_direction = HomogeneusCoor.normalize(
                            Vec3d.positions_to_vec(pixel_pos_in_space, self.camera.eye))
                        ray = Ray(start_pos=pixel_pos_in_space, direction=ray_direction)
                        color += self.calculate_color_rt(200, ray, first_ray=1)
                #
                color = color * (1 / self.subsample)

            pixel_list[y % num_thread][x] = color.fix_range().getRGB()

        return pixel_list


    def calculate_color_rt(self, nof_bounce, coming_ray, max_bounce=NOF_MAX_BOUNCE, inside_ray=False, first_ray = 0 ):
        color = Color(0, 0, 0, 0)

        if nof_bounce == 0 or max_bounce == 0:
            return color
        else:

            closest_node, intersect, c2, i2 = self.__find_closest_obj(coming_ray)

            if closest_node != -1:
                point = intersect.p0
                surface_normal = closest_node.normal_of_a_point(point)
                ray_dir = coming_ray.direction
                bias = 0.00001

                for light in self.scene.lights:
                    if light.light_type == "point":
                        if self.__is_that_closest_to_light(closest_node, point, light):
                            lightDir = Vec3d.positions_to_vec(light.position.normalize(), point.normalize()).normalize()
                            nDotL = max(surface_normal.normalize().dot_product(lightDir), 0.0)
                            spec = 0

                            #Calculate Specular
                            if first_ray:
                                viewDir = (ray_dir*(-1)).normalize()
                                halfwayDir = lightDir.add(viewDir).normalize()
                                spec = pow(max(halfwayDir.dot_product(surface_normal), 0), 1/0.01)

                            frasnel_value = 1 if closest_node.material.refractive == 0 else Ray.calculate_frasnel(
                                ray_dir, surface_normal, object_ior=closest_node.material.ior)
                            reflective_rate = closest_node.material.reflective

                            # REFRACTIVE CALCULATION
                            if closest_node.material.refractive != 0:
                                point= point.add(surface_normal* (-bias))
                                ray = Ray.calculated_bounced_ray(ray_dir, surface_normal, point)
                                reflective_color_value1 = self.calculate_color_rt(nof_bounce - 1, ray, max_bounce - 1)

                                refract_ray_tmp = Ray.calculate_refract_ray(ray_dir, surface_normal, closest_node.material.ior,point)

                                p1 = closest_node.intersect(refract_ray_tmp).p1
                                tmp_normal = closest_node.normal_of_a_point(p1)

                                refract_ray = Ray.calculate_refract_ray(refract_ray_tmp.direction, tmp_normal * (-1), 1 / closest_node.material.ior, p1)
                                reflect_ray = Ray.calculated_bounced_ray(refract_ray_tmp.direction, tmp_normal * (-1), p1)

                                refractive_color_value = Color(0, 0, 0, 0)
                                reflective_color_value = Color(0, 0, 0, 0)
                                if inside_ray == True:
                                    return Color(0, 0, 0)
                                elif refract_ray.direction != 0:
                                    refractive_color_value = self.calculate_color_rt(nof_bounce - 1, refract_ray, max_bounce - 1, inside_ray= True)
                                    reflective_color_value = self.calculate_color_rt(nof_bounce - 1, reflect_ray, max_bounce = max_bounce - 1)

                                else:
                                    reflective_color_value = self.calculate_color_rt(nof_bounce - 1, refract_ray_tmp, max_bounce - 1)

                                color += ((((refractive_color_value) * (1 - frasnel_value) + reflective_color_value * (
                                 frasnel_value)) * (1 - reflective_rate) + reflective_color_value1 * (reflective_rate)
                                 )) * (1 - frasnel_value +spec ) * (light.intensity)

                            # REFLECTIVE CALCULATION
                            elif reflective_rate != 0:
                                point = point.add(surface_normal * (bias))
                                ray = Ray.calculated_bounced_ray(ray_dir, surface_normal, point)
                                reflective_color_value = self.calculate_color_rt(nof_bounce - 1, ray, max_bounce - 1)

                                color += (closest_node.material.color * (1 -reflective_rate) + reflective_color_value * (
                                        reflective_rate)) * (nDotL + spec) * light.intensity

                            else:
                                color += closest_node.material.color * (1 - reflective_rate) * (nDotL + spec) * light.intensity

                    elif light.light_type == "dome":
                        new_color = self.calculate_dome_light(intersect.p0, closest_node, light)
                        color += Color(new_color[0], new_color[1], new_color[2]) * light.intensity

            return color



    def calculate_dome_light(self, intersect_point, closest_node, dome_light):
        nof_intersect = 0
        point = intersect_point
        sample = self.ambient_sample
        for r in range(sample):
            sample_ray = closest_node.generate_random_hemisphere_ray(point)

            for n in self.scene.nodes:
                if n != closest_node:
                    int_sample = n.intersect(sample_ray)
                    if (int_sample != -1):
                        if (int_sample.t0 > 0):
                            nof_intersect += 1
                            break
        if nof_intersect == 0:
            new_color = dome_light.color.getRGB()

        else:
            new_color = [x * ((sample - nof_intersect) / sample) for x in
                         closest_node.material.color.getRGB()]

        '''
        light_value = 0
        point = intersect_point  # ustte tanimli sil sonra
        sample = self.ambient_sample
        for r in range(sample):
            sample_ray = closest_node.generate_random_hemisphere_ray(point)

            tmp_nodes = copy.deepcopy(self.scene.nodes)
            tmp_nodes.append(dome_light.lightObj)

            tmp_closest = self.__find_closest_obj(sample_ray, tmp_nodes)
            if tmp_closest == dome_light.lightObj:
                light_value += 1
        if light_value == 0:
            new_color = dome_light.color.getRGB()

        else:
            new_color = [x * ((light_value) / sample) for x in
                         closest_node.material.color.getRGB()]
        '''
        return new_color


    def reflective_bounce(self, nof_refra_bounce, nof_refle_bounce, nof_bounce, coming_ray, surface_normal,
                          intersect_pos):
        return 0


    def __find_closest_obj(self, ray, nodes = None):
        if nodes == None:
            nodes = self.scene.nodes

        shortest = -1
        intersect = -1
        closest_node = -1

        closest_node_2 = -1
        intersect_2 = -1
        for node in nodes:

            intersect_tmp = node.intersect(ray)
            if (intersect_tmp != -1):
                if (intersect_tmp.t0 >= 0):
                    d = Pos3d.cal_distance_value(intersect_tmp.p0, ray.start_pos)
                    if (shortest == -1) or (shortest > d):
                        closest_node_2 = closest_node
                        intersect_2 = intersect
                        shortest = d
                        closest_node = node
                        intersect = intersect_tmp

        return closest_node, intersect, closest_node_2, intersect_2


    def __is_that_closest_to_light(self, object, intersect_point, light):
        ray_from_light = Ray(start_pos=light.position,
                             direction=HomogeneusCoor.normalize(
                                 Vec3d.positions_to_vec(intersect_point, light.position)))
        closest_node_to_light, tmp_i, c2, i2 = self.__find_closest_obj(ray_from_light)
        if object == closest_node_to_light:
            if tmp_i.p0.compare(intersect_point):
                return True
            if closest_node_to_light.material.refractive > 0:
                return True

        elif closest_node_to_light.material.refractive > 0:
            if object == c2:
                if i2.p0.compare(intersect_point):
                    return True
        return False


    def __is_that_closest_to_light_with_check_light(self, object, intersect_point, light):
        ray_from_light = Ray(start_pos=light.position,
                             direction=HomogeneusCoor.normalize(
                                 Vec3d.positions_to_vec(intersect_point, light.position)))
        closest_node_to_light, tmp_i, c2, i2 = self.__find_closest_obj(ray_from_light)
        if object == closest_node_to_light:
            if tmp_i.p0.compare(intersect_point):
                return True, False
            if closest_node_to_light.material.refractive > 0:
                return True, True

        elif closest_node_to_light.material.refractive > 0:
            if object == c2:
                if i2.p0.compare(intersect_point):
                    return True, True
        return False, False
