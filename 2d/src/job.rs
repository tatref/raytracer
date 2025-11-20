use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub enum Message {
    WorkRequest { world_id: usize, columns: Vec<i32> },
}
