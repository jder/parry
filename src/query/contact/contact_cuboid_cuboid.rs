use crate::math::{Pose, Real};
#[cfg(feature = "dim2")]
use crate::query::details::clip_segment_segment_with_normal;
#[cfg(feature = "dim3")]
use crate::query::PointQuery;
use crate::query::{sat, Contact};
use crate::shape::Cuboid;

/// Contact between two cuboids.
#[inline]
pub fn contact_cuboid_cuboid(
    pos12: &Pose,
    cuboid1: &Cuboid,
    cuboid2: &Cuboid,
    prediction: Real,
) -> Option<Contact> {
    let pos21 = pos12.inverse();

    let sep1 = sat::cuboid_cuboid_find_local_separating_normal_oneway(cuboid1, cuboid2, pos12);
    if sep1.0 > prediction {
        return None;
    }

    let sep2 = sat::cuboid_cuboid_find_local_separating_normal_oneway(cuboid2, cuboid1, &pos21);
    if sep2.0 > prediction {
        return None;
    }

    #[cfg(feature = "dim3")]
    let sep3 = sat::cuboid_cuboid_find_local_separating_edge_twoway(cuboid1, cuboid2, pos12);
    #[cfg(feature = "dim3")]
    if sep3.0 > prediction {
        return None;
    }

    #[cfg(feature = "dim2")]
    {
        let best_sep = if sep2.0 > sep1.0 {
            (sep2.0, pos12.rotation * -sep2.1)
        } else {
            sep1
        };
        let normal2 = pos12.rotation.inverse() * -best_sep.1;
        let feature1 = cuboid1.support_feature(best_sep.1);
        let feature2 = cuboid2.support_feature(normal2);
        let (clip_a, clip_b) = clip_segment_segment_with_normal(
            (feature1.vertices[0], feature1.vertices[1]),
            (pos12 * feature2.vertices[0], pos12 * feature2.vertices[1]),
            best_sep.1,
        )?;
        let dist_a = (clip_a.1 - clip_a.0).dot(best_sep.1);
        let dist_b = (clip_b.1 - clip_b.0).dot(best_sep.1);
        let (point1, point2_world, dist) = if dist_a <= dist_b {
            (clip_a.0, clip_a.1, dist_a)
        } else {
            (clip_b.0, clip_b.1, dist_b)
        };

        return Some(Contact::new(
            point1,
            pos12.inverse_transform_point(point2_world),
            best_sep.1,
            normal2,
            dist,
        ));
    }

    #[cfg(feature = "dim3")]
    {
        // The best separating axis is face-vertex.
        if sep1.0 >= sep2.0 && sep1.0 >= sep3.0 {
            let pt2_1 = cuboid2.support_point(pos12, -sep1.1);
            let proj1 = cuboid1.project_local_point(pt2_1, false);

            let separation = (pt2_1 - proj1.point).dot(sep1.1);
            let (mut normal1, mut dist) = (pt2_1 - proj1.point).normalize_and_length();

            if separation < 0.0 || dist <= Real::default_epsilon() {
                normal1 = sep1.1;
                dist = separation;
            }

            if dist > prediction {
                return None;
            }

            return Some(Contact::new(
                proj1.point,
                pos12.inverse_transform_point(pt2_1),
                normal1,
                pos12.rotation.inverse() * -normal1,
                dist,
            ));
        }

        if sep2.0 >= sep1.0 && sep2.0 >= sep3.0 {
            let pt1_2 = cuboid1.support_point(&pos21, -sep2.1);
            let proj2 = cuboid2.project_local_point(pt1_2, false);

            let separation = (pt1_2 - proj2.point).dot(sep2.1);
            let (mut normal2, mut dist) = (pt1_2 - proj2.point).normalize_and_length();

            if separation < 0.0 || dist <= Real::default_epsilon() {
                normal2 = sep2.1;
                dist = separation;
            }

            if dist > prediction {
                return None;
            }

            return Some(Contact::new(
                pos12 * pt1_2,
                proj2.point,
                pos12.rotation * -normal2,
                normal2,
                dist,
            ));
        }

        use crate::query::{details, ClosestPoints};
        let edge1 = cuboid1.local_support_edge_segment(sep3.1);
        let edge2 = cuboid2.local_support_edge_segment(pos21.rotation * -sep3.1);

        match details::closest_points_segment_segment(pos12, &edge1, &edge2, prediction) {
            ClosestPoints::Disjoint => None,
            ClosestPoints::WithinMargin(a, b) => {
                let normal1 = sep3.1;
                let normal2 = pos12.rotation.inverse() * -normal1;
                Some(Contact::new(a, b, normal1, normal2, sep3.0))
            }
            ClosestPoints::Intersecting => unreachable!(),
        }
    }
}
