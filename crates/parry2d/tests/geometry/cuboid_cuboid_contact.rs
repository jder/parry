use parry2d::{
    math::{Pose, Rotation, Vector},
    query::{self, ClosestPoints},
    shape::Cuboid,
};

#[test]
fn cuboid_cuboid_contact_handles_half_turn_negative_zero() {
    let cuboid1 = Cuboid::new(Vector::new(190.0, 110.0));
    let cuboid2 = Cuboid::new(Vector::new(21.0, 12.0));

    let pos1 = Pose::IDENTITY;
    let pos2 = Pose::from_parts(
        Vector::new(206.69629, -115.64551),
        Rotation::from_cos_sin_unchecked(-1.0, -0.0),
    );

    assert!(query::intersection_test(&pos1, &cuboid1, &pos2, &cuboid2).unwrap());
    assert_eq!(
        query::closest_points(&pos1, &cuboid1, &pos2, &cuboid2, 10.0).unwrap(),
        ClosestPoints::Intersecting
    );
    assert!(
        query::contact(&pos1, &cuboid1, &pos2, &cuboid2, 10.0)
            .unwrap()
            .is_some()
    );
}
