use parry2d::{
    math::{Pose, Rotation, Vector},
    query::{self, ClosestPoints},
    shape::Cuboid,
};

#[test]
fn cuboid_cuboid_contact_matches_sampled_closest_points_queries() {
    let cuboid1 = Cuboid::new(Vector::new(1.0, 1.0));
    let cuboid2 = Cuboid::new(Vector::new(1.0, 1.0));
    let pos1 = Pose::IDENTITY;
    let prediction = 0.35;

    for &angle in &[0.0, 0.1, 0.2, 0.4, 0.7, 1.0] {
        for xi in 10..=30 {
            for yi in 10..=30 {
                let x = xi as f32 * 0.1;
                let y = yi as f32 * 0.1;
                let pos2 = Pose::from_parts(Vector::new(x, y), Rotation::new(angle));
                let closest =
                    query::closest_points(&pos1, &cuboid1, &pos2, &cuboid2, prediction).unwrap();
                let contact = query::contact(&pos1, &cuboid1, &pos2, &cuboid2, prediction)
                    .unwrap();

                match (closest, contact) {
                    (ClosestPoints::Disjoint, Some(contact)) => panic!(
                        "disjoint but contact: angle={angle}, x={x}, y={y}, dist={}",
                        contact.dist
                    ),
                    (ClosestPoints::WithinMargin(_, _), None) => panic!(
                        "within-margin but no contact: angle={angle}, x={x}, y={y}"
                    ),
                    (_, Some(contact)) if contact.dist > prediction => panic!(
                        "contact beyond prediction: angle={angle}, x={x}, y={y}, dist={}",
                        contact.dist
                    ),
                    _ => {}
                }
            }
        }
    }
}
