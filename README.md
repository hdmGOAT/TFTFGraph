# TFTF Edge: Revolutionizing Jeepney Routes & Fare Optimization 🚗💨

## 🚀 Overview

Welcome to **TFTF Edge**, a state-of-the-art C++ application designed to optimize jeepney routing and fare calculations. TFTF Edge uses a **Temporal Flexible Transfer and Fare Graph (TFTF)** to empower commuters and jeepney operators alike with smart, real-time route and fare recommendations. 

No more waiting in traffic or paying unpredictable fares—TFTF Edge makes jeepney commuting **faster**, **smarter**, and **affordable**!

---

## 🌟 Key Features

- **Dynamic Route Modeling**: Route flexibility at its best! TFTF Edge calculates the most efficient jeepney routes without strict schedules, adapting in real-time to traffic and other variables.
- **Fare Optimization**: Never overpay! Fare calculations take into account distance, transfer points, and jeepney density to give you the most cost-effective options.
- **Efficient Transfers**: Find the best transfer points automatically, ensuring that your journey is quick and smooth.
- **Density-Aware Routing**: Routes consider jeepney density, helping you avoid crowded rides for a more comfortable commute.
- **Real-Time Adjustments**: Get real-time route and fare updates based on live traffic and jeepney conditions.

---

## ⚡ How It Works

### **1. Route Calculation**
TFTF Edge dynamically computes optimal jeepney routes between two stations, factoring in transfers, traffic conditions, and jeepney density.

### **2. Fare Calculation**
Once the optimal route is determined, TFTF Edge calculates the fare based on the route’s distance, number of transfers, and other dynamic conditions.

### **3. Dynamic Adjustments**
The app continuously updates routes and fares, ensuring you're always getting the best option in real-time.

**Under the hood**, we use a **TFTF Graph**: 
- **Nodes** represent jeepney stations.
- **Edges** are the possible routes and transfers between these stations.
- The weight of the edges dynamically adjusts based on **time**, **traffic**, and **jeepney density**.

---

## 🧠 The Math Behind the Magic

Here's how we calculate the most optimal routes and fares:

### **Route Cost Calculation**

The formula to select routes looks like this:

```
Total Cost = Transfer Cost × Density Factor
```

This calculates how efficient the route is in terms of transfer points and how crowded the jeepney is.

### **Fare Calculation**

The fare is based on:

```
TFare = (BASEFARE × N) + (DISTANCE / 1000 × FPKM)
```

Where:
- `BASEFARE`: Base fare of the jeepney ride.
- `N`: The number of transfers.
- `DISTANCE`: The total distance covered.
- `FPKM`: Fare per kilometer.

This model provides accurate fare estimates and helps with financial planning for your trip.

---

## 🛠 Installation

### Prerequisites:
- **C++17** or later
- **CMake** for building the project
- **Linux** or **Windows** development environment

### Steps to Get Started:

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/yourusername/tftf-edge.git
   ```

2. **Navigate to the Project Folder**:

   ```bash
   cd tftf-edge
   ```

3. **Build the Project** using CMake:

   ```bash
   mkdir build
   cd build
   cmake ..
   make
   ```

4. **Run the Application**:

   ```bash
   ./tftf-edge
   ```

Now you can start exploring dynamic jeepney routing and fare calculation in action!

---

## 🤝 Contributing

**TFTF Edge** is a community-driven project, and we welcome contributions from anyone interested in improving the world of jeepney routing and fare optimization! Here's how you can help:

1. **Fork** the repository and clone it to your local machine.
2. **Make your changes**—whether it’s adding a feature or fixing a bug.
3. **Submit a Pull Request** with a description of what you've done.

---

## 📜 License

TFTF Edge is open-source and licensed under the **MIT License**. See the [LICENSE](LICENSE) file for more details.

---

## 🎉 Let's Build the Future of Public Transport Together!

**TFTF Edge** is more than just an app—it’s a movement towards smarter, more flexible jeepney routes and fairer fare systems. Whether you're a commuter, developer, or transport operator, you’re helping to shape the future of public transportation.
