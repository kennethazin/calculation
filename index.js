Cesium.Ion.defaultAccessToken = import.meta.env.VITE_CESIUM_ACCESS_TOKEN;

const viewer = new Cesium.Viewer("cesiumContainer", {
  terrain: Cesium.Terrain.fromWorldTerrain(),
  shouldAnimate: false,
});

const startTime = Cesium.JulianDate.fromIso8601("1969-07-16T13:32:00Z");
const satelliteStopTime = Cesium.JulianDate.fromIso8601("1969-07-16T18:05:20Z");

async function initialize() {
  const dataSource = await viewer.dataSources.add(
    Cesium.CzmlDataSource.load("./saturn_v_trajectory_orientation.czml")
  );

  const satellite = dataSource.entities.getById("SaturnV");
  if (satellite) {
    satellite.viewFrom = new Cesium.Cartesian3(-300, 20, 100);
  } else {
    console.error(
      "Satellite entity with ID '3304422' not found in the CZML file."
    );
  }

  // Set the camera to follow the satellite by default
  viewer.trackedEntity = satellite;

  // Particle system for thrusters
  const thrusterParticles = new Cesium.ParticleSystem({
    image: "./thruster_particle.png", // Path to particle image
    startColor: Cesium.Color.RED.withAlpha(0.7),
    endColor: Cesium.Color.YELLOW.withAlpha(0.3),
    startScale: 1.0,
    endScale: 4.0,
    minimumParticleLife: 0.5,
    maximumParticleLife: 1.5,
    minimumSpeed: 5.0,
    maximumSpeed: 10.0,
    emissionRate: 50,
    emitter: new Cesium.ConeEmitter(Cesium.Math.toRadians(30)),
    modelMatrix: Cesium.Matrix4.IDENTITY,
    lifetime: 16.0,
  });

  viewer.scene.primitives.add(thrusterParticles);

  function updateThrusterParticles() {
    const currentTime = Cesium.JulianDate.toDate(
      viewer.clock.currentTime
    ).getTime();

    // Define burn times (in milliseconds since epoch)
    const burnTimes = [
      {
        start: Date.parse("1969-07-16T13:32:00Z"),
        end: Date.parse("1969-07-16T13:34:48Z"),
      }, // Stage 1 (168 seconds)
      {
        start: Date.parse("1969-07-16T13:34:48Z"),
        end: Date.parse("1969-07-16T13:40:54Z"),
      }, // Stage 2 (366 seconds)
      {
        start: Date.parse("1969-07-16T13:40:54Z"),
        end: Date.parse("1969-07-16T13:43:18Z"),
      }, // Stage 3 Burn 1 (144 seconds)
      {
        start: Date.parse("1969-07-16T16:59:18Z"),
        end: Date.parse("1969-07-16T17:05:54Z"),
      }, // Stage 3 Burn 2 (336 seconds)
    ];

    // Check if current time is within any burn period
    const isBurning = burnTimes.some(
      (burn) => currentTime >= burn.start && currentTime <= burn.end
    );

    thrusterParticles.show = isBurning;

    if (satellite && isBurning) {
      // Ensure position and orientation exist before accessing them
      const position = satellite.position?.getValue(viewer.clock.currentTime);
      const orientation = satellite.orientation?.getValue(
        viewer.clock.currentTime
      );

      if (position && orientation) {
        const modelMatrix = Cesium.Transforms.headingPitchRollToFixedFrame(
          position,
          Cesium.HeadingPitchRoll.fromQuaternion(orientation)
        );
        thrusterParticles.modelMatrix = modelMatrix;
      }
    }
  }

  viewer.clock.onTick.addEventListener(updateThrusterParticles);

  const viewModel = {
    show: true,
    intensity: 2.0,
    distortion: 10.0,
    dispersion: 0.4,
    haloWidth: 0.4,
    dirtAmount: 0.4,
  };

  const lensFlare = viewer.scene.postProcessStages.add(
    Cesium.PostProcessStageLibrary.createLensFlareStage()
  );

  function updatePostProcess() {
    lensFlare.enabled = Boolean(viewModel.show);
    lensFlare.uniforms.intensity = Number(viewModel.intensity);
    lensFlare.uniforms.distortion = Number(viewModel.distortion);
    lensFlare.uniforms.ghostDispersal = Number(viewModel.dispersion);
    lensFlare.uniforms.haloWidth = Number(viewModel.haloWidth);
    lensFlare.uniforms.dirtAmount = Number(viewModel.dirtAmount);
    lensFlare.uniforms.earthRadius = Cesium.Ellipsoid.WGS84.maximumRadius;

    // Increase the resolution of the lens flare reflection
    lensFlare.uniforms.resolution = 1024; // Set a higher resolution value
  }
  updatePostProcess();

  // Add event listeners for buttons
  document.getElementById("satelliteButton").addEventListener("click", () => {
    viewer.clock.stopTime = satelliteStopTime;
    viewer.clock.currentTime = startTime;
    viewer.clock.multiplier = 30;
    viewer.timeline.zoomTo(startTime, satelliteStopTime);
    viewer.trackedEntity = satellite; // Ensure the camera follows the satellite
  });

  document.getElementById("trackingAuto").addEventListener("click", () => {
    satellite.trackingReferenceFrame = Cesium.TrackingReferenceFrame.AUTODETECT;
  });

  document.getElementById("trackingInertial").addEventListener("click", () => {
    satellite.trackingReferenceFrame = Cesium.TrackingReferenceFrame.INERTIAL;
  });

  document.getElementById("trackingVelocity").addEventListener("click", () => {
    satellite.trackingReferenceFrame = Cesium.TrackingReferenceFrame.VELOCITY;
  });

  document.getElementById("trackingENU").addEventListener("click", () => {
    satellite.trackingReferenceFrame = Cesium.TrackingReferenceFrame.ENU;
  });
}

initialize();
