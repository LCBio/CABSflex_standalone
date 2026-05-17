// This script runs early (not deferred) to set scroll restoration
if ('scrollRestoration' in history) {
  history.scrollRestoration = 'manual';
}

// Ensure we are at the top as early as possible
window.scrollTo(0, 0);

// Reset scroll position to top before leaving the page
// This ensures that if the browser tries to "restore" the scroll position 
// on the next page load (even in manual mode), it restores "0".
window.addEventListener('beforeunload', () => {
  window.scrollTo(0, 0);
});

document.addEventListener('DOMContentLoaded', () => {
  const toggle = document.querySelector(".nav-toggle");
  const sidebar = document.getElementById("sidebar");

  if (toggle && sidebar) {
    toggle.addEventListener("click", () => {
      const open = sidebar.classList.toggle("is-open");
      toggle.setAttribute("aria-expanded", open ? "true" : "false");
    });
  }

  // Handle links to top or internal anchors
  document.addEventListener('click', (e) => {
    const link = e.target.closest('a');
    if (!link) return;
    const href = link.getAttribute('href');
    if (!href) return;

    if (href === '#' || href === '#top' || href === '#site-top') {
      e.preventDefault();
      window.scrollTo({ top: 0, left: 0, behavior: 'smooth' });
      history.pushState(null, null, ' ');
    }
  });

  // Final fallback for page load
  window.addEventListener('load', () => {
    if (!window.location.hash || window.location.hash === '#' || window.location.hash === '#top') {
      window.scrollTo(0, 0);
    }
  });

  // --- Target Highlight Visual Pulse ---
  const highlightTarget = () => {
    const hash = window.location.hash;
    if (hash) {
      try {
        const decodedHash = decodeURIComponent(hash);
        const target = document.querySelector(decodedHash);
        if (target) {
          target.classList.remove('target-highlight');
          void target.offsetWidth; // Trigger reflow to restart animation
          target.classList.add('target-highlight');
          
          // Also highlight next paragraph if it's a reference page section
          if (target.tagName === 'H2' && target.nextElementSibling && target.nextElementSibling.tagName === 'P') {
            const nextPara = target.nextElementSibling;
            nextPara.classList.remove('target-highlight');
            void nextPara.offsetWidth;
            nextPara.classList.add('target-highlight');
          }
        }
      } catch (err) {
        console.error('Error in anchor highlighting:', err);
      }
    }
  };

  window.addEventListener('hashchange', highlightTarget);
  // Run on load with a slight delay to ensure layout is complete and scroll position is settled
  if (window.location.hash) {
    setTimeout(highlightTarget, 300);
  }

  // --- Reference Hover Previews ---
  let referencesDoc = null;
  const fetchReferences = async () => {
    if (referencesDoc) return referencesDoc;
    try {
      // Find relative path to references.html (handles different pages in root docs dir)
      const response = await fetch('references.html');
      if (!response.ok) throw new Error('Failed to fetch references.html');
      const htmlText = await response.text();
      const parser = new DOMParser();
      referencesDoc = parser.parseFromString(htmlText, 'text/html');
      return referencesDoc;
    } catch (err) {
      console.error('Error loading reference previews:', err);
      return null;
    }
  };

  // Create tooltip element
  const tooltip = document.createElement('div');
  tooltip.className = 'reference-tooltip';
  document.body.appendChild(tooltip);

  let activeTimeout = null;

  document.addEventListener('mouseover', async (e) => {
    const link = e.target.closest('a');
    if (!link) return;
    const href = link.getAttribute('href');
    if (!href || !href.includes('references.html#')) return;

    const hashIndex = href.indexOf('#');
    const hash = decodeURIComponent(href.substring(hashIndex));
    const doc = await fetchReferences();
    if (!doc) return;

    const targetHeading = doc.querySelector(hash);
    if (!targetHeading) return;

    // Extract citation title and body paragraph
    const headingText = targetHeading.textContent.replace('¶', '').trim();
    const nextPara = targetHeading.nextElementSibling;
    const paraHTML = nextPara && nextPara.tagName === 'P' ? nextPara.innerHTML : '';

    tooltip.innerHTML = `<strong>${headingText}</strong><hr style="margin:8px 0; border:0; border-top:1px solid rgba(15, 98, 254, 0.15);">${paraHTML}`;
    
    // Position tooltip relative to hovered link
    const rect = link.getBoundingClientRect();
    const tooltipWidth = 340;
    tooltip.style.width = `${tooltipWidth}px`;
    
    // Pre-show with opacity 0 to calculate offsetHeight
    tooltip.style.left = '-9999px';
    tooltip.style.top = '-9999px';
    tooltip.classList.add('is-visible');
    const tooltipHeight = tooltip.offsetHeight;
    tooltip.classList.remove('is-visible');

    let left = rect.left + window.scrollX + (rect.width / 2) - (tooltipWidth / 2);
    let top = rect.top + window.scrollY - tooltipHeight - 12;

    // Check bounds & flip tooltip below the link if it overflows the top viewport boundary
    if (left < 12) left = 12;
    if (left + tooltipWidth > window.innerWidth - 12) {
      left = window.innerWidth - tooltipWidth - 12;
    }
    if (rect.top - tooltipHeight - 12 < 12) {
      top = rect.bottom + window.scrollY + 12;
      tooltip.classList.add('is-below');
    } else {
      tooltip.classList.remove('is-below');
    }

    tooltip.style.left = `${left}px`;
    tooltip.style.top = `${top}px`;

    clearTimeout(activeTimeout);
    tooltip.classList.add('is-visible');
  });

  document.addEventListener('mouseout', (e) => {
    const link = e.target.closest('a');
    if (!link) return;
    const href = link.getAttribute('href');
    if (!href || !href.includes('references.html#')) return;

    activeTimeout = setTimeout(() => {
      tooltip.classList.remove('is-visible');
    }, 200);
  });
});
